#include <TApplication.h>
#include <TCanvas.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <TRandom3.h>
#include <memory>
#include <chrono>
#include <numeric>
#include <cmath>
#include <thread>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMicroscopicTypes.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/RandomEngine.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Random.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/ComponentGrid.hh"
// FFTW and threading for Poisson solver
#include <fftw3.h>
#include <mutex>
#include <future>


using namespace Garfield;

// ----------------------------- PoissonFFT2D class ------------------------------
// 2D Poisson solver using FFTW. Accepts a 2D charge density (x,z)
// and returns Ex and Ez across the same Nx x Nz grid. The solver
// uses zero-padding to reduce wrap-around and multi-threading
// for per-k operations.
class PoissonFFT2D {
public:
  PoissonFFT2D(int Nx_, int Nz_, double spacing_, int pad = 10, int nthreads = 0)
      : Nx(Nx_), Nz(Nz_), spacing(spacing_), Npad(pad) {
    // padded sizes
    Nx_pad = Nx + 2 * Npad;
    Nz_pad = Nz + 2 * Npad;
    Nz_r2c = Nz_pad / 2 + 1;
    if (nthreads > 0) {
      nthreads_ = nthreads;
    } else {
      unsigned int hc = std::thread::hardware_concurrency();
      nthreads_ = (hc > 0) ? static_cast<int>(hc) : 1;
    }

    // allocate arrays
    rho_in = (double*)fftw_malloc(sizeof(double) * Nx_pad * Nz_pad);
    rho_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx_pad * Nz_r2c);
    phi_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx_pad * Nz_r2c);
    Ex_fft  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx_pad * Nz_r2c);
    Ez_fft  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx_pad * Nz_r2c);
    Ex_real = (double*)fftw_malloc(sizeof(double) * Nx_pad * Nz_pad);
    Ez_real = (double*)fftw_malloc(sizeof(double) * Nx_pad * Nz_pad);

    // create plans
    forward_plan = fftw_plan_dft_r2c_2d(Nx_pad, Nz_pad, rho_in, rho_fft, FFTW_MEASURE);
    backward_plan_Ex = fftw_plan_dft_c2r_2d(Nx_pad, Nz_pad, Ex_fft, Ex_real, FFTW_MEASURE);
    backward_plan_Ez = fftw_plan_dft_c2r_2d(Nx_pad, Nz_pad, Ez_fft, Ez_real, FFTW_MEASURE);
  }
  ~PoissonFFT2D() {
    if (forward_plan) fftw_destroy_plan(forward_plan);
    if (backward_plan_Ex) fftw_destroy_plan(backward_plan_Ex);
    if (backward_plan_Ez) fftw_destroy_plan(backward_plan_Ez);
    if (rho_in) fftw_free(rho_in);
    if (rho_fft) fftw_free(rho_fft);
    if (phi_fft) fftw_free(phi_fft);
    if (Ex_fft) fftw_free(Ex_fft);
    if (Ez_fft) fftw_free(Ez_fft);
    if (Ex_real) fftw_free(Ex_real);
    if (Ez_real) fftw_free(Ez_real);
  }

  // Solve Poisson and fill Ex and Ez (both Nx*Nz flattened vectors)
  void Solve(const std::vector<double>& rho, std::vector<double>& Ex_out, std::vector<double>& Ez_out) {
    if ((int)rho.size() != Nx * Nz) {
      std::cerr << "PoissonFFT2D::Solve: rho size mismatch\n";
      return;
    }
    // zero input
    const int total_pad = Nx_pad * Nz_pad;
    std::fill(rho_in, rho_in + total_pad, 0.0);

    // copy to padded input (centered)
    for (int ix = 0; ix < Nx; ++ix) {
      for (int iz = 0; iz < Nz; ++iz) {
        int ixp = ix + Npad;
        int izp = iz + Npad;
        rho_in[ixp * Nz_pad + izp] = rho[ix * Nz + iz];
      }
    }

    // forward transform
    fftw_execute(forward_plan);

    // compute phi_fft = rho_fft / (eps0 * k2) in parallel and Ex/Ez in k domain
    const double Lx = Nx_pad * spacing;
    const double Lz = Nz_pad * spacing;

    auto work = [&](int ix_start, int ix_end) {
      for (int ix = ix_start; ix < ix_end; ++ix) {
        int kx = (ix <= Nx_pad / 2) ? ix : ix - Nx_pad; // kx = ix for the pad, else its the actual value we want
        double kx_val = 2.0 * Pi * kx / Lx;
        for (int iz = 0; iz < Nz_r2c; ++iz) {
          int kz = iz; // r2c output limited to positive side
          double kz_val = 2.0 * Pi * kz / Lz;
          int idx = ix * Nz_r2c + iz;
          const double a = rho_fft[idx][0];
          const double b = rho_fft[idx][1];
          double k2 = kx_val * kx_val + kz_val * kz_val;
          if (k2 < 1e-20) { // sets small values of k to 0 to minimise noise
            phi_fft[idx][0] = 0.0; phi_fft[idx][1] = 0.0;
            Ex_fft[idx][0] = 0.0; Ex_fft[idx][1] = 0.0;
            Ez_fft[idx][0] = 0.0; Ez_fft[idx][1] = 0.0;
          } else {
            // phi = rho/(eps0*k2), complex division (rho is complex)
            double denom = VacuumPermittivity * k2;
            // phi_fft = (a + i b) / denom
            double phi_re = a / denom;
            double phi_im = b / denom;
            phi_fft[idx][0] = phi_re;
            phi_fft[idx][1] = phi_im;

            // Ex_fft = -i kx * phi = -i kx*(phi_re + i phi_im) = kx*phi_im + i*(-kx*phi_re)
            Ex_fft[idx][0] = kx_val * phi_im;
            Ex_fft[idx][1] = -kx_val * phi_re;

            // Ez_fft = -i kz * phi
            Ez_fft[idx][0] = kz_val * phi_im;
            Ez_fft[idx][1] = -kz_val * phi_re;
          }
        }
      }
    };

    // Launch threads
    std::vector<std::thread> threads;
    int per_thread = Nx_pad / nthreads_;  // the number of x values each thread does
    int start = 0;
    for (int t = 0; t < nthreads_; ++t) {
      int end = (t == nthreads_ - 1) ? Nx_pad : start + per_thread; // the last thread goes to the end
      threads.emplace_back(work, start, end);
      start = end;
    }
    for (auto &th : threads) th.join();

    // inverse transforms back to real domain
    fftw_execute(backward_plan_Ex);
    fftw_execute(backward_plan_Ez);

    // normalize and copy back central region
    const double Ntotal = (double)(Nx_pad * Nz_pad);
    Ex_out.assign(Nx * Nz, 0.0);
    Ez_out.assign(Nx * Nz, 0.0);
    for (int ix = 0; ix < Nx; ++ix) {
      for (int iz = 0; iz < Nz; ++iz) {
        int ixp = ix + Npad;
        int izp = iz + Npad;
        int idx_pad = ixp * Nz_pad + izp;  // flattened array indexing
        // note backward produced scaled by Ntotal, so divide
        Ex_out[ix * Nz + iz] = Ex_real[idx_pad] / Ntotal;
        Ez_out[ix * Nz + iz] = Ez_real[idx_pad] / Ntotal;
      }
    }
  }

private:
  int Nx, Nz;
  int Nx_pad, Nz_pad, Nz_r2c, Npad;
  double spacing;
  int nthreads_;
  double *rho_in; 
  fftw_complex *rho_fft, *phi_fft, *Ex_fft, *Ez_fft;
  double *Ex_real, *Ez_real;
  fftw_plan forward_plan;
  fftw_plan backward_plan_Ex;
  fftw_plan backward_plan_Ez;
};
// ----------------------------- End PoissonFFT2D --------------------------------


bool readTransferFunction(Sensor& sensor) {
  std::ifstream infile;
  infile.open("mdt_elx_delta.txt", std::ios::in);
  if (!infile) {
    std::cerr << "Could not read delta response function.\n";
    return false;
  }
  std::vector<double> times;
  std::vector<double> values;
  while (!infile.eof()) {
    double t = 0., f = 0.;
    infile >> t >> f;
    if (infile.eof() || infile.fail()) break;
    times.push_back(1.e3 * t);
    values.push_back(f);
  }
  infile.close();
  sensor.SetTransferFunction(times, values);
  return true;
}
/*
double firstPassageTime(double t0, double D, double a) {
  double u = RndmUniform();
  double z = std::erfinv(1-u);
  return t0 + (a*a)/(4 * D * z*z);
}
*/
int main(int argc, char* argv[]) {

auto time_start = std::chrono::steady_clock::now();
// Read seed before creating TApplication
int seed = 123456;
if (argc > 1) {
  seed = std::stoi(argv[1]);
  std::cout << "Using user-specified seed: " << seed << std::endl;
} else {
  std::cout << "Using default seed: " << seed << std::endl;
}

// Set ROOT seed
gRandom->SetSeed(seed);

// create RNG engine and explicitly seed it
Garfield::RandomEngineRoot randomEngine; //default constructor
randomEngine.SetSeed(seed);
Garfield::Random::SetEngine(randomEngine);

  TApplication app("app", &argc, argv);

  // Make a gas medium.
  MediumMagboltz gas;
  gas.LoadGasFile("n2_78.08_o2_20.95_ar_0.93_co2_0.04_1atm.gas");
  gas.LoadIonMobility("IonMobility_N2+_N2.txt");
  gas.LoadNegativeIonMobility("NegIonMobility_O2-_air.txt");
 
  // Make a component with analytic electric field.
  ComponentAnalyticField cmp;

  const double dt = 100.; //time between loops (dt > tstep)
  const double tstep = 100.; //monte-carlo step size (ns)
  const double v_drift = 220e-9; //cm/ns
  const double tmin = -0.5 * tstep; 
  const std::size_t nbins = 3000;
  const bool stop_at_max_time = true;
  const size_t max_time = nbins*tstep;

  // Grid parameters.
  const int Nx = 100;
  const int Ny = 1;
  const int Nz = 100;
  const double spacing = v_drift * tstep * 10; //(cm)
  //const double spacing = 0.0001;
  const double xgrid = Nx * spacing;
  const double zgrid = Nz * spacing;
  const double alpha = 1.72e-15; // recombination coefficient (cm^3/ns)
  const bool RecordRecombinationPositions = false;

  // Plate Separation (cm)
  const double xMin = -0.15, xMax = 0.15;
  // Width in beam direction
  const double yMin = -0.15, yMax = 0.15;
  // Length? Of the plates (vertically)
  const double zMin = -0.15, zMax = 0.15;
  
  const double vAnode = 15.;
  const double vCathode = 0.;
  
 cmp.AddPlaneX(xMax, vAnode, "detector"); // plane that hosts our detector
 cmp.AddPlaneX(xMin, vCathode);

 //cmp.AddStripOnPlaneX('z', xMax, yMin, yMax, "detector");
 cmp.SetMedium(&gas);

  ComponentGrid grid;
  grid.SetMesh(Nx, Ny, Nz, -xgrid/2, xgrid/2, yMin,
                 yMax, -zgrid/2, zgrid/2);
  grid.SetUniformElectricField(0., 0., 0.); 
  grid.SetMedium(&gas);

  // Make a sensor.
  Sensor sensor(&cmp);
  sensor.AddElectrode(&cmp, "detector");
  sensor.AddComponent(&grid);

  // option to make strip sensors for positional resolution
  const bool stripSensor = false;
  const int nStrips = 20;
  const double stripWidth = (yMax - yMin) / nStrips;
  sensor.SetTimeWindow(tmin, tstep, nbins);
  std::vector<std::string> stripNames;
  if (stripSensor) {
    for (double y_0 = yMin; y_0 <= (yMin + (nStrips - 1) * stripWidth); y_0 += stripWidth ) {
      stripNames.push_back(std::to_string(y_0));
      cmp.AddStripOnPlaneX('z', xMax, y_0, y_0 + stripWidth, stripNames.back());
      sensor.AddElectrode(&cmp, stripNames.back());
    }
  }
  // Set the delta reponse function.
  if (!readTransferFunction(sensor)) return 0;

  bool JumpIonsOutsideGridToPlate = false;
  if (JumpIonsOutsideGridToPlate) {sensor.SetArea(-xgrid/2, yMin, -zgrid/2, xgrid/2, yMax, zgrid/2);}
  else {sensor.SetArea(xMin, yMin, zMin, xMax, yMax, zMax);}
    
  // Set up Heed.
  TrackHeed track(&sensor);
  track.SetParticle("proton");
  track.SetEnergy(250.e6 + ProtonMass);


  AvalancheMC drift;
    drift.SetSensor(&sensor);
    drift.EnableSignalCalculation();
    drift.SetTimeSteps(tstep);
    drift.EnableDensityMap();
    drift.EnableRecombination(true, alpha);
    const bool SpaceCharge = true;

  AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);


  TCanvas* cD = nullptr;
  ViewDrift driftView;
  constexpr bool plotDrift = false;  //can cause massive memory gain over time
  if (plotDrift) {
    cD = new TCanvas("cD", "", 600, 600);
    driftView.SetCanvas(cD);
    drift.EnablePlotting(&driftView);
    track.EnablePlotting(&driftView);
  }

  TCanvas* cS = nullptr;
  constexpr bool plotSignal = true;
  if (plotSignal) {
    cS = new TCanvas("cS", "", 600, 600);
    }  // Two GUI Windows

  const double x0 = 0;
  const double y0 = yMin;
  const double z0 = 0;
  const std::size_t nTracks = 500;
  const int multiplicity = 10;  // charges per ion/electron
  double recombine_num = 0;
  std::vector<std::vector<double>> ion_recombination_positions;
  std::vector<std::vector<double>> negion_recombination_positions;
  std::vector<AvalancheMC::EndPoint> IonsLeftGrid;
  std::vector<AvalancheMC::EndPoint> NegativeIonsLeftGrid;
  int ElectronsLeftGridCount = 0;
  
  
  
  sensor.ClearSignal();
  for (std::size_t j = 0; j < nTracks; ++j) {
    double offset = randomEngine.Draw();
    double x_proton = x0 + (offset - 0.5) * 0.000; //spread over 0.01 mm
    track.NewTrack(x_proton, y0, z0, 0, 0, 1, 0);
    for (const auto& cluster : track.GetClusters()) {
      //remove clusters that are unphysically out of the detector
      if (cluster.y < yMin || cluster.y > yMax) {
        continue;
      }
       //adds the particles to our drifting functions 
      for (const auto& Ion : cluster.ions) {
        drift.AddIon(Ion.x, Ion.y, Ion.z, Ion.t, multiplicity);
      }
      for (const auto& electron : cluster.electrons) {
        aval.AddElectron(electron.x, electron.y, electron.z, electron.t, 0., 0., 0., 0., multiplicity); 
      }
    }
  }
    std::cout <<"Initial Electrons: " << aval.GetElectrons().size() << std::endl;
    std::cout <<"Initial Ions: " << drift.GetIons().size() << std::endl;

  //for (double t = 0; t < (nbins * tstep); t += dt) { original time based loop
  double t = 0;
  size_t particleNum = aval.GetElectrons().size() + drift.GetIons().size() + drift.GetNegativeIons().size();
  
  // Pre-allocate arrays and constants used by SpaceCharge computations
  std::vector<double> chargeDensity;
  std::function<int(int,int)> idx;
  double xGridMin = -xgrid/2 + spacing/2; // centre of the first cell in x/z
  double zGridMin = -zgrid/2 + spacing/2;
  std::unique_ptr<PoissonFFT2D> poissonSolver;
  if (SpaceCharge) {
    chargeDensity.resize(Nx * Nz);
    idx = [Nz](int ix, int iz) { return ix * Nz + iz; };
    int Npad = 10;   // padding on each side for the FFT poisson solver
    // initialize Poisson solver, 0 threads auto allocates
    poissonSolver = std::make_unique<PoissonFFT2D>(Nx, Nz, spacing, /*pad=*/Npad, /*threads=*/0);
  }

  
  // ----------------------------- main while loop -----------------------------------
  while (particleNum > 0) {
    if (stop_at_max_time && t > max_time) {break;}

    // handle electron drift and attachment
    if (!aval.GetElectrons().empty()) {
      aval.SetTimeWindow(t, t + dt);
      aval.ResumeAvalanche(); //drift the electrons only if there are some left
      // check for electron attachment and add negative ions.
      for (const auto& electron : aval.GetElectrons()) {
        if (electron.status == -7) {
            const auto& p1 = electron.path.back();
            drift.AddNegativeIon(p1.x, p1.y, p1.z, p1.t, multiplicity);
        }
        if (electron.status == StatusLeftDriftArea && JumpIonsOutsideGridToPlate) {
            ElectronsLeftGridCount += 1;
        }
      }
    }
    
    grid.ClearFields();  // clear old densities/fields
    if (t > 0 && SpaceCharge) {
      bool loaded = grid.LoadElectricField("electric_field.xyz", "XYZ", false, false, 1.0, 1.0, 1.0);
    } else {
      grid.SetUniformElectricField(0., 0., 0.);
    }

    // add positive ions to the grid TODO stop the particles outside the grid running because we can check faster
    for (const auto& ion : drift.GetIons()) {
      if (!ion.path.empty()) {
        const auto& p1 = ion.path.back();
        grid.AddIon(p1.x, p1.y, p1.z, multiplicity);
      }
    }

    // add negative ions to the grid
    for (const auto& negion : drift.GetNegativeIons()) {
      if (!negion.path.empty()) {
        const auto& p1 = negion.path.back();
        grid.AddNegativeIon(p1.x, p1.y, p1.z, multiplicity);
      }
    }

    // -------------------- Space Charge -------------------------
    if (SpaceCharge) {
      // reset pre-allocated arrays and reuse precomputed constants
      std::fill(chargeDensity.begin(), chargeDensity.end(), 0.0);
      
      // make the charge density map (flattened array)
      const double y = 0.0;
      for (int ix = 0; ix < Nx; ++ix) {
        double x = xGridMin + ix * spacing;
        for (int iz = 0; iz < Nz; ++iz) {
          double z = zGridMin + iz * spacing;
          double rhoIon = 0.0, rhoNegIon = 0.0;
          grid.IonDensity(x, y, z, rhoIon);
          grid.NegativeIonDensity(x, y, z, rhoNegIon);

          // Net charge density per voxel
          chargeDensity[idx(ix, iz)] = ElementaryCharge * (rhoIon - rhoNegIon); // fC/cm^3
        }
      }

      // get the electric field from the charge density
      std::vector<double> Ex(Nx * Nz, 0.0);
      std::vector<double> Ez(Nx * Nz, 0.0);
      if (poissonSolver) {
        poissonSolver->Solve(chargeDensity, Ex, Ez);
      } else {
        // fall back to 0 field
      }

      // Save Ex and Ez to a file in XYZ format for ComponentGrid to read
      std::ofstream efieldFile("electric_field.xyz");
      if (!efieldFile) {
        std::cerr << "Cannot open file for writing electric field.\n";
      } else {
        efieldFile << std::scientific << std::setprecision(6);
        const double y = 0.0; // single slice in y
        for (int ix = 0; ix < Nx; ++ix) {
          double x = xGridMin + ix * spacing;
          for (int iz = 0; iz < Nz; ++iz) {
              double z = zGridMin + iz * spacing;
              double ex = Ex[idx(ix, iz)];
              double ey = 0.0; // assume no Ey
              double ez = Ez[idx(ix, iz)];
              efieldFile << x << " " << y << " " << z << " "
                        << ex << " " << ey << " " << ez << "\n";
          }
        }
        efieldFile.close();
        std::cout << "Electric field saved to electric_field.xyz\n";
      } 
    }
// ------------------------- End Space Charge ------------------------------------
    drift.SetTimeWindow(t, t + dt);
    // drift the positive and negative ions
    drift.ResumeAvalanche();

    std::cout << "Negative Ions: " << drift.GetNegativeIons().size() << std::endl;
    std::cout << "Positive Ions: " << drift.GetIons().size() << std::endl;
    // record recombined particles
    for (auto &ion : drift.GetIons()) {
                if (ion.status == -9) {
                  const auto& p1 = ion.path.back();
                  ion_recombination_positions.push_back({p1.x, p1.y, p1.z});
                  recombine_num += 1;
                }
                else if (ion.status == StatusLeftDriftArea && JumpIonsOutsideGridToPlate) {
                  IonsLeftGrid.push_back(ion);
                }
            }
    for (auto &negion : drift.GetNegativeIons()) {
                if (negion.status == -9) {
                  const auto& p1 = negion.path.back();
                  negion_recombination_positions.push_back({p1.x, p1.y, p1.z});
                  recombine_num += 1;
                }
                else if (negion.status == StatusLeftDriftArea && JumpIonsOutsideGridToPlate) {
                  NegativeIonsLeftGrid.push_back(negion);
                }
            }

    std::cout << t + dt << "ns simulated" << std::endl;
    t += dt;
    particleNum = aval.GetElectrons().size() + drift.GetIons().size() + drift.GetNegativeIons().size();
  }


  std::cout << "Recombined particles : " << recombine_num << std::endl;
  if (ElectronsLeftGridCount > 0) {std::cout << "WARNING " << ElectronsLeftGridCount << " electrons left the grid" << std::endl;}  
  
if (plotDrift) {
  cD->Clear();
  cmp.PlotCell(cD);
  constexpr bool twod = true;
  constexpr bool drawaxis = false;
  driftView.Plot(twod, drawaxis);
}
  

//sensor.ConvoluteSignals();

// option to display integrated signal
bool integrateSignal = true;
if (integrateSignal) {
  sensor.IntegrateSignal("detector");
  double total_charge = sensor.GetSignal("detector", nbins - 1);
  std::cout << "Total collected charge: " << total_charge << " fC" << std::endl;
  std::cout << "Corresponding number of electrons: " << total_charge / ElementaryCharge << std::endl;
}

// below here is all outputting data to files and plotting
  bool saveSignal = false;
  if (saveSignal) { std::ofstream outfile;
    outfile.open("signal.txt", std::ios::out);
    for (unsigned int i = 0; i < nbins; ++i) {
    const double t = (i + 0.5) * tstep;
    const double f = sensor.GetSignal("detector", i);
    const double fe = sensor.GetElectronSignal("detector", i);
    const double fh = sensor.GetIonSignal("detector", i);
    outfile << t << " " << f << " " << fe << " " << fh << "\n";
    }
    outfile.close();
  }

  if (stripSensor) {
    std::ofstream outfile;
    sensor.IntegrateSignals();
    outfile.open("strip_signals.txt", std::ios::out);
    outfile << "Strip_y_position(cm) total_charge(fC)\n";
    for (const auto& name : stripNames) {
      const double strip_total_charge = sensor.GetSignal(name, nbins - 1);
      outfile << name << " " << strip_total_charge <<"\n";
    }
  }

  if (RecordRecombinationPositions) {
    std::ofstream outfile;
    outfile.open("recombination_positions.txt", std::ios::out);
    outfile << "Ion Recombination Positions (cm):\n";
    for (const auto& pos : ion_recombination_positions) {
      outfile << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }
    outfile.close();
    outfile.open("negion_recombination_positions.txt", std::ios::out);
    outfile << "Negative Ion Recombination Positions (cm):\n";
    for (const auto& pos : negion_recombination_positions) {
      outfile << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }
    outfile.close();
  }
  
  if (plotSignal) { 
    sensor.PlotSignal("detector", cS);
  }

  bool saveIonPositions = false;
  if (saveIonPositions) {
    std::ofstream outfile;
    std::vector<std::vector<double>> ion_positions;
    outfile.open("ion_positions.txt", std::ios::out);
    outfile << "x(cm) y(cm) z(cm) t=" << t / 1000 << "\u03BCs\n";
    for (const auto& ion : drift.GetIons()) {
      const auto& p1 = ion.path.back();
      ion_positions.push_back({p1.x, p1.y, p1.z});
      outfile << p1.x << " " << p1.y << " " << p1.z << "\n";
    }
    outfile.close();
    outfile.open("negion_positions.txt", std::ios::out);
    outfile << "x(cm) y(cm) z(cm) t=" << t / 1000 << "\u03BCs\n";
    for (const auto& negion : drift.GetNegativeIons()) {
      const auto& p1 = negion.path.back();
      ion_positions.push_back({p1.x, p1.y, p1.z});
      outfile << p1.x << " " << p1.y << " " << p1.z << "\n";
    }
    outfile.close();
  }

  // timer
  auto time_end = std::chrono::steady_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
  std::cout << "Total Signal: " << sensor.GetSignal("detector", nbins - 1) << " Ion Signal: " << sensor.GetIonSignal("detector", nbins-1) << " Electron + Negion Signal: " << sensor.GetElectronSignal("detector", nbins-1) << std::endl;
  std::cout << "Elapsed time: " << elapsed_ms << " ms" << std::endl;


  app.Run(kTRUE);
}