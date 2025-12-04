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

class PoissonFFT3D {
public:
  PoissonFFT3D(int Nx_, int Ny_, int Nz_, double sx_, double sy_, double sz_, int pad = 10, int nthreads = 0)
      : Nx(Nx_), Ny(Ny_), Nz(Nz_), sx(sx_), sy(sy_), sz(sz_), Npad(pad) {
    Nx_pad = Nx + 2 * Npad;
    Ny_pad = Ny + 2 * Npad;
    Nz_pad = Nz + 2 * Npad;
    Nz_r2c = Nz_pad / 2 + 1;
    if (nthreads > 0) nthreads_ = nthreads;
    else { unsigned int hc = std::thread::hardware_concurrency(); nthreads_ = (hc > 0) ? static_cast<int>(hc) : 1; }

    const size_t in_size = static_cast<size_t>(Nx_pad) * Ny_pad * Nz_pad;
    const size_t fftc_size = static_cast<size_t>(Nx_pad) * Ny_pad * Nz_r2c;
    rho_in = (double*)fftw_malloc(sizeof(double) * in_size);
    rho_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftc_size);
    phi_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftc_size);
    Ex_fft  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftc_size);
    Ey_fft  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftc_size);
    Ez_fft  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftc_size);
    Ex_real = (double*)fftw_malloc(sizeof(double) * in_size);
    Ey_real = (double*)fftw_malloc(sizeof(double) * in_size);
    Ez_real = (double*)fftw_malloc(sizeof(double) * in_size);

    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads_);
    forward_plan = fftw_plan_dft_r2c_3d(Nx_pad, Ny_pad, Nz_pad, rho_in, rho_fft, FFTW_MEASURE);
    backward_plan_Ex = fftw_plan_dft_c2r_3d(Nx_pad, Ny_pad, Nz_pad, Ex_fft, Ex_real, FFTW_MEASURE);
    backward_plan_Ey = fftw_plan_dft_c2r_3d(Nx_pad, Ny_pad, Nz_pad, Ey_fft, Ey_real, FFTW_MEASURE);
    backward_plan_Ez = fftw_plan_dft_c2r_3d(Nx_pad, Ny_pad, Nz_pad, Ez_fft, Ez_real, FFTW_MEASURE);
  }
  ~PoissonFFT3D() {
    if (forward_plan) fftw_destroy_plan(forward_plan);
    if (backward_plan_Ex) fftw_destroy_plan(backward_plan_Ex);
    if (backward_plan_Ey) fftw_destroy_plan(backward_plan_Ey);
    if (backward_plan_Ez) fftw_destroy_plan(backward_plan_Ez);
    if (rho_in) fftw_free(rho_in);
    if (rho_fft) fftw_free(rho_fft);
    if (phi_fft) fftw_free(phi_fft);
    if (Ex_fft) fftw_free(Ex_fft);
    if (Ey_fft) fftw_free(Ey_fft);
    if (Ez_fft) fftw_free(Ez_fft);
    if (Ex_real) fftw_free(Ex_real);
    if (Ey_real) fftw_free(Ey_real);
    if (Ez_real) fftw_free(Ez_real);
    fftw_cleanup_threads();
  }

  void Solve(const std::vector<double>& rho, std::vector<double>& Ex_out, std::vector<double>& Ey_out, std::vector<double>& Ez_out) {
    if ((int)rho.size() != Nx * Ny * Nz) {
      std::cerr << "PoissonFFT3D::Solve: rho size mismatch\n";
      return;
    }
    const size_t total_pad = static_cast<size_t>(Nx_pad) * Ny_pad * Nz_pad;
    std::fill(rho_in, rho_in + total_pad, 0.0);
    for (int ix = 0; ix < Nx; ++ix) {
      for (int iy = 0; iy < Ny; ++iy) {
        for (int iz = 0; iz < Nz; ++iz) {
          int ixp = ix + Npad;
          int iyp = iy + Npad;
          int izp = iz + Npad;
          size_t idx_pad = static_cast<size_t>(ixp) * Ny_pad * Nz_pad + static_cast<size_t>(iyp) * Nz_pad + izp;
          size_t idx_src = static_cast<size_t>(ix) * Ny * Nz + static_cast<size_t>(iy) * Nz + iz;
          rho_in[idx_pad] = rho[idx_src];
        }
      }
    }

    fftw_execute(forward_plan);
    const double eps0 = Garfield::VacuumPermittivity;
    const double Lx = Nx_pad * sx;
    const double Ly = Ny_pad * sy;
    const double Lz = Nz_pad * sz;
    auto work = [&](int ix_start, int ix_end) {
      for (int ix = ix_start; ix < ix_end; ++ix) {
        int kx = (ix <= Nx_pad / 2) ? ix : ix - Nx_pad;
        double kx_val = 2.0 * Garfield::Pi * kx / Lx;
        for (int iy = 0; iy < Ny_pad; ++iy) {
          int ky = (iy <= Ny_pad / 2) ? iy : iy - Ny_pad;
          double ky_val = 2.0 * Garfield::Pi * ky / Ly;
          for (int iz = 0; iz < Nz_r2c; ++iz) {
            int kz = iz;
            double kz_val = 2.0 * Garfield::Pi * kz / Lz;
            size_t idx = static_cast<size_t>(ix) * Ny_pad * Nz_r2c + static_cast<size_t>(iy) * Nz_r2c + iz;
            const double a = rho_fft[idx][0];
            const double b = rho_fft[idx][1];
            double k2 = kx_val * kx_val + ky_val * ky_val + kz_val * kz_val;
            if (k2 < 1e-20) {
              phi_fft[idx][0] = 0.0; phi_fft[idx][1] = 0.0;
              Ex_fft[idx][0] = 0.0; Ex_fft[idx][1] = 0.0;
              Ey_fft[idx][0] = 0.0; Ey_fft[idx][1] = 0.0;
              Ez_fft[idx][0] = 0.0; Ez_fft[idx][1] = 0.0;
            } else {
              double denom = eps0 * k2;
              double phi_re = a / denom;
              double phi_im = b / denom;
              phi_fft[idx][0] = phi_re;
              phi_fft[idx][1] = phi_im;
              Ex_fft[idx][0] = kx_val * phi_im;
              Ex_fft[idx][1] = -kx_val * phi_re;
              Ey_fft[idx][0] = ky_val * phi_im;
              Ey_fft[idx][1] = -ky_val * phi_re;
              Ez_fft[idx][0] = kz_val * phi_im;
              Ez_fft[idx][1] = -kz_val * phi_re;
            }
          }
        }
      }
    };

    std::vector<std::thread> threads;
    int per_thread = Nx_pad / nthreads_;
    int start = 0;
    for (int t = 0; t < nthreads_; ++t) {
      int end = (t == nthreads_ - 1) ? Nx_pad : start + per_thread;
      threads.emplace_back(work, start, end);
      start = end;
    }
    for (auto &th : threads) th.join();

    fftw_execute(backward_plan_Ex);
    fftw_execute(backward_plan_Ey);
    fftw_execute(backward_plan_Ez);

    const double Ntotal = static_cast<double>(Nx_pad) * Ny_pad * Nz_pad;
    Ex_out.assign(Nx * Ny * Nz, 0.0);
    Ey_out.assign(Nx * Ny * Nz, 0.0);
    Ez_out.assign(Nx * Ny * Nz, 0.0);
    for (int ix = 0; ix < Nx; ++ix) {
      for (int iy = 0; iy < Ny; ++iy) {
        for (int iz = 0; iz < Nz; ++iz) {
          int ixp = ix + Npad;
          int iyp = iy + Npad;
          int izp = iz + Npad;
          size_t idx_pad = static_cast<size_t>(ixp) * Ny_pad * Nz_pad + static_cast<size_t>(iyp) * Nz_pad + izp;
          size_t idx_out = static_cast<size_t>(ix) * Ny * Nz + static_cast<size_t>(iy) * Nz + iz;
          Ex_out[idx_out] = Ex_real[idx_pad] / Ntotal;
          Ey_out[idx_out] = Ey_real[idx_pad] / Ntotal;
          Ez_out[idx_out] = Ez_real[idx_pad] / Ntotal;
        }
      }
    }
  }

private:
  int Nx, Ny, Nz;
  int Nx_pad, Ny_pad, Nz_pad, Nz_r2c, Npad;
  double sx, sy, sz;
  int nthreads_;
  double *rho_in; 
  fftw_complex *rho_fft, *phi_fft, *Ex_fft, *Ey_fft, *Ez_fft;
  double *Ex_real, *Ey_real, *Ez_real;
  fftw_plan forward_plan;
  fftw_plan backward_plan_Ex;
  fftw_plan backward_plan_Ey;
  fftw_plan backward_plan_Ez;
};
// ----------------------------- End PoissonFFT3D --------------------------------

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
  const double tmin = -0.5 * tstep; 
  const std::size_t nbins = 100;
  const bool stop_at_max_time = true;
  const size_t max_time = nbins*tstep;
 //-----------------------------------------Geometry-----------------------------------------------
 // proton travels along the y axis.
 // Plate Separation (cm)
  const double xMin = -15., xMax = 15.;
  // Width in beam direction
  const double yMin = -0.25, yMax = 0.25;
  // Length? Of the plates (vertically)
  const double zMin = -20., zMax = 20.;
  
  const double vAnode = 0.;
  const double vCathode = -1500.;
  const double E_mag = std::abs(vCathode/(yMax-yMin));

  cmp.SetMedium(&gas);
  // cmp.AddPixelOnPlaneY(yMax ,-13.125 ,13.125 , -17.5, 17.5, "detector"); // define size of detector makes the program run really slowly???
  cmp.AddPlaneY(yMax, vAnode, "detector"); // plane that hosts our detector
  cmp.AddPlaneY(yMin, vCathode); // supposed to have two high voltage planes but analytic field only allows 2.
  // Grid parameters.
  double vx_negion,vy_negion,vz_negion;
  cmp.GetMedium(0,0,0)->NegativeIonVelocity(0,E_mag,0,0,0,0,vx_negion,vy_negion,vz_negion);
  const double v_drift = 8.12176e-06; //cm/ns O2- drift velocity in air at 3000 V/cm
  const double guide_spacingy = -vy_negion * tstep * 5; //(cm) We define the spacing as 10 average drift lengths.
  const double spacing_transverse = 0.001; //cm
  const int Nx = 20; //number of grid spaces in y
  const int Ny = std::round((yMax-yMin) / guide_spacingy);
  const int Nz = 20;
  const double spacingy = (yMax - yMin)/ Ny;
  const double xgrid = Nx * spacing_transverse;
  const double zgrid = Nz * spacing_transverse;
  const double alpha = 1.72e-15; // recombination coefficient (cm^3/ns)
  const bool RecordRecombinationPositions = false;
  std::cout << "Grid spans: x=+-" << xgrid/2 << " z=+-" << zgrid/2 << " and has Ny=" << Ny << std::endl;
  
// mesh
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
  // NEEDS UPDATING FOR NEW DETECTOR GEOMETRY TODO
  const bool stripSensor = false;
  const int nStrips = 20;
  const double stripWidth = (yMax - yMin) / nStrips;
  std::vector<std::string> stripNames;
  if (stripSensor) {
    for (double y_0 = yMin; y_0 <= (yMin + (nStrips - 1) * stripWidth); y_0 += stripWidth ) {
      stripNames.push_back(std::to_string(y_0));
      cmp.AddStripOnPlaneY('z', xMax, y_0, y_0 + stripWidth, stripNames.back());
      sensor.AddElectrode(&cmp, stripNames.back());
    }
  }
  sensor.SetTimeWindow(tmin, tstep, nbins);
  
  sensor.SetArea(xMin, yMin, zMin, xMax, yMax, zMax); //particles that leave the area are removed from simulation
    
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
  constexpr bool plotDrift = true;  //can cause massive memory gain over time
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

  const double x0 = 0; // centre of start of proton beam in x
  const double sigmax = 0.46;
  const double sigmaz = 0.66;
  const double y0 = yMin;
  const double z0 = 0; // centre of start of proton beam in z
  const std::size_t nTracks = 500; // number of protons
  const int multiplicity = 1;  // charges per ion/electron
  double recombine_num = 0; // number of recombinations so far (counter)
  std::vector<std::vector<double>> ion_recombination_positions;
  std::vector<std::vector<double>> negion_recombination_positions;
  std::vector<AvalancheMC::EndPoint> IonsLeftGrid;
  std::vector<AvalancheMC::EndPoint> NegativeIonsLeftGrid;
  int ElectronsLeftGridCount = 0;
  
  
  
  sensor.ClearSignal();
  for (std::size_t j = 0; j < nTracks; ++j) {
    double x_proton = RndmGaussian(x0, sigmax);
    double z_proton = RndmGaussian(z0, sigmaz);
    track.NewTrack(x_proton, y0, z_proton, 0, 0, 1, 0);
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
  
  // Pre-allocate arrays and constants used by SpaceCharge computations (3D)
  std::vector<double> chargeDensity;
  std::function<int(int,int,int)> idx3d;
  double xGridMin = -xgrid/2 + spacing_transverse/2; // centre of the first cell in x and z
  double zGridMin = -zgrid/2 + spacing_transverse/2;
  double yGridMin = yMin + spacingy/2;  // centre of first cell in y
  std::unique_ptr<PoissonFFT3D> poissonSolver;
  if (SpaceCharge) {
    chargeDensity.resize(Nx * Ny * Nz);
    idx3d = [Ny, Nz](int ix, int iy, int iz) { return (ix * Ny + iy) * Nz + iz; }; // indexes the 1D arrays as 3D.
    int Npad = 10;   // padding on each side for the FFT poisson solver
    // initialize 3D Poisson solver, 0 threads auto allocates
    poissonSolver = std::make_unique<PoissonFFT3D>(Nx, Ny, Nz, spacing_transverse, spacingy, spacing_transverse, /*pad=*/Npad, /*threads=*/0);
  }

  
  // ----------------------------- main while loop -----------------------------------
  while (particleNum > 0) {
    if (stop_at_max_time && t >= max_time) {break;}

    // handle electron drift and attachment
    if (!aval.GetElectrons().empty()) {
      aval.SetTimeWindow(t, t + dt);
      aval.ResumeAvalanche();
      std::cout << "Yep!" << std::endl; //drift the electrons only if there are some left
      // check for electron attachment and add negative ions.
      for (const auto& electron : aval.GetElectrons()) {
        if (electron.status == -7) {
            const auto& p1 = electron.path.back();
            drift.AddNegativeIon(p1.x, p1.y, p1.z, p1.t, multiplicity);
        }
        if (electron.status == StatusLeftDriftArea) {
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

    // -------------------- Space Charge -----------------------------------------
    if (SpaceCharge) {
      // reset pre-allocated arrays and reuse precomputed constants
      std::fill(chargeDensity.begin(), chargeDensity.end(), 0.0);
      
      // make the charge density map (flattened array) over x,y,z
      size_t n_ions = drift.GetIons().size();
      size_t n_negions = drift.GetNegativeIons().size();
      for (int ix = 0; ix < Nx; ++ix) {
        double x = xGridMin + ix * spacing_transverse;
        for (int iy = 0; iy < Ny; ++iy) {
          double y = yGridMin + iy * spacingy;
          for (int iz = 0; iz < Nz; ++iz) {
            double z = zGridMin + iz * spacing_transverse;
            double rhoIon = 0.0, rhoNegIon = 0.0;
            if (n_ions > 0) {
              grid.IonDensity(x, y, z, rhoIon);
            }
            if (n_negions > 0) {
              grid.NegativeIonDensity(x, y, z, rhoNegIon);
            }
            if (rhoIon == false) {rhoIon = 0;}
            if (rhoNegIon == false) {rhoNegIon = 0;}
            // Net charge density per voxel
            chargeDensity[idx3d(ix, iy, iz)] = Garfield::ElementaryCharge * (rhoIon - rhoNegIon); // fC/cm^3
          }
        }
      }

      // get the electric field from the charge density
      std::vector<double> Ex(Nx * Ny * Nz, 0.0);
      std::vector<double> Ey(Nx * Ny * Nz, 0.0);
      std::vector<double> Ez(Nx * Ny * Nz, 0.0);
      if (poissonSolver) {
        bool hasCharge = std::any_of(chargeDensity.begin(), chargeDensity.end(), [](double v){ return std::abs(v) > 1e-30; });
        if (hasCharge) poissonSolver->Solve(chargeDensity, Ex, Ey, Ez);
      } else {
        // fall back to 0 field
      }

      // Save Ex and Ez to a file in XYZ format for ComponentGrid to read
      std::ofstream efieldFile("electric_field.xyz");
      if (!efieldFile) {
        std::cerr << "Cannot open file for writing electric field.\n";
      } else {
        efieldFile << std::scientific << std::setprecision(6);
        for (int ix = 0; ix < Nx; ++ix) {
          double x = xGridMin + ix * spacing_transverse;
          for (int iy = 0; iy < Ny; ++iy) {
            double y = yGridMin + iy * spacingy;
            for (int iz = 0; iz < Nz; ++iz) {
              double z = zGridMin + iz * spacing_transverse;
              double ex = Ex[idx3d(ix, iy, iz)];
              double ey = Ey[idx3d(ix, iy, iz)];
              double ez = Ez[idx3d(ix, iy, iz)];
              efieldFile << x << " " << y << " " << z << " "
                        << ex << " " << ey << " " << ez << "\n";
            }
          }
        }
        efieldFile.close();
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
                else if (ion.status == StatusLeftDriftArea) {
                  IonsLeftGrid.push_back(ion);
                }
            }
    for (auto &negion : drift.GetNegativeIons()) {
                if (negion.status == -9) {
                  const auto& p1 = negion.path.back();
                  negion_recombination_positions.push_back({p1.x, p1.y, p1.z});
                  recombine_num += 1;
                }
                else if (negion.status == StatusLeftDriftArea) {
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
  

// option to display integrated signal
bool integrateSignal = true;
if (integrateSignal) {
  sensor.IntegrateSignal("detector");
  double total_charge = sensor.GetSignal("detector", nbins - 1);
  std::cout << "Total collected charge: " << total_charge << " fC" << std::endl;
  std::cout << "Corresponding number of electrons: " << total_charge / Garfield::ElementaryCharge << std::endl;
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