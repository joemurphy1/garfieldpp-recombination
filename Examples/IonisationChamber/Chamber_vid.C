#include <TApplication.h>
#include <TCanvas.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <chrono>
#include <numeric>
#include <filesystem>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/RandomEngine.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Random.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/ComponentGrid.hh"


using namespace Garfield;


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
 

  const double dt = 100.; //time between loops (dt > tstep)
  const double tstep = 100.; //monte-carlo step size (ns)
  const double v_drift = 220e-9; //cm/ns

    // Plate Seperation (cm)
  const double xMin = -0.15, xMax = 0.15;
  // Width in beam direction
  const double yMin = -0.15, yMax = 0.15;
  // Length? Of the plates (vertically)
  const double zMin = -0.15, zMax = 0.15;
  // Voltages
  const double vAnode = 15.;
  const double vCathode = 0.;

  // Grid parameters.
  const int Nx = 100;
  const int Ny = 1;
  const int Nz = 100;
  const double spacing = v_drift * dt * 10; //(cm)
  const double xgrid = Nx * spacing;
  const double zgrid = Nz * spacing;
  const double alpha = 1.72e-15; // recombination coefficient (cm^3/ns)
  const bool RecordRecombinationPositions = false;


  ComponentGrid grid;
  grid.SetMedium(&gas);
  grid.SetMesh(Nx, Ny, Nz, -xgrid/2, xgrid/2, yMin,
                 yMax, -zgrid/2, zgrid/2);
  grid.SetUniformElectricField(0., 0., 0.); 
  grid.SetMedium(&gas);

  // Make a component with analytic electric field.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
 cmp.AddPlaneX(xMin, vCathode);
 cmp.AddPlaneX(xMax, vAnode, "anode");
 cmp.AddStripOnPlaneX('z', xMax, yMin, yMax, "detector");
  // Make a sensor.
  Sensor sensor(&cmp);
  sensor.AddElectrode(&cmp, "detector");
  sensor.AddComponent(&grid);

  // option to make strip sensors for positional resolution
  const bool stripSensor = false;
  const int nStrips = 20;
  const double stripWidth = (yMax - yMin) / nStrips;
  std::vector<std::string> stripNames;
  if (stripSensor) {
    for (double y_0 = yMin; y_0 <= (yMin + (nStrips - 1) * stripWidth); y_0 += stripWidth ) {
      stripNames.push_back(std::to_string(y_0));
      cmp.AddStripOnPlaneX('z', xMax, y_0, y_0 + stripWidth, stripNames.back());
      sensor.AddElectrode(&cmp, stripNames.back());
    }
  }

  // Set the signal time window.
  const double tmin = -0.5 * tstep; 
  const std::size_t nbins = 10000;
  const bool stop_at_max_time = true;
  const size_t max_time = nbins*tstep;
  sensor.SetTimeWindow(tmin, tstep, nbins);
  // Set the delta reponse function.
  if (!readTransferFunction(sensor)) return 0;
  sensor.ClearSignal();

  bool generateVideo = true;
  const int videoInterval = 10; //interval in number of dt steps
  if (generateVideo) {
    std::filesystem::path parent = "particle_positions";

    try {
        // Create the parent directory
        std::filesystem::create_directories(parent);

        // Create the subdirectories
        std::filesystem::create_directories(parent / "electrons");
        std::filesystem::create_directories(parent / "ions");
        std::filesystem::create_directories(parent / "negions");

        std::cout << "Directories created successfully.\n";

    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << "\n";
    }
  }

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
  const std::size_t nTracks = 1000;
  double recombine_num = 0;
  std::vector<std::vector<double>> ion_recombination_positions;
  std::vector<std::vector<double>> negion_recombination_positions;
  
  
  
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
        drift.AddIon(Ion.x, Ion.y, Ion.z, Ion.t);
      }
      for (const auto& electron : cluster.electrons) {
        aval.AddElectron(electron.x, electron.y, electron.z, electron.t, 0., 0., 0., 0.); 
      }
    }
  }
    std::cout <<"Initial Electrons: " << aval.GetElectrons().size() << std::endl;
    std::cout <<"Initial Ions: " << drift.GetIons().size() << std::endl;

  //for (double t = 0; t < (nbins * tstep); t += dt) { original time based loop
  double t = 0;
  size_t particleNum = aval.GetElectrons().size() + drift.GetIons().size() + drift.GetNegativeIons().size();
  
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
              drift.AddNegativeIon(p1.x, p1.y, p1.z, p1.t, 1);
          }
        }
      }
    
    grid.ClearFields();  // clear old densities/fields
    if (t > 0) {
      bool loaded = grid.LoadElectricField("electric_field.xyz", "XYZ", false, false, 1.0, 1.0, 1.0);
    } else {
      grid.SetUniformElectricField(0., 0., 0.);
    }

    // add positive ions to the grid
    const int multiplicity = 1;  // or however many charges per ion you want
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
    std::vector<std::vector<double>> chargeDensity(Nx, std::vector<double>(Nz, 0.0));
    double y = 0.;
    for (int ix = 0; ix < Nx; ++ix) {
      double x = -xgrid/2 + spacing/2 + ix * spacing;
      for (int iz = 0; iz < Nz; ++iz) {
        double z = -zgrid/2 + spacing/2 + iz * spacing;
        double rhoIon = 0.0, rhoNegIon = 0.0;
        grid.IonDensity(x, y, z, rhoIon);
        grid.NegativeIonDensity(x, y, z, rhoNegIon);

        // Net charge density per voxel
        chargeDensity[ix][iz] = ElementaryCharge * (rhoIon - rhoNegIon); // fC/cm^3
      }
    }

    std::vector<std::vector<double>> Ex(Nx, std::vector<double>(Nz, 0.0));
    std::vector<std::vector<double>> Ez(Nx, std::vector<double>(Nz, 0.0));
    const double k = 1.0 / (FourPiEpsilon0); // [cm·fC^-1] in CGS-like units if rho in fC/cm^3

    for (int ix = 0; ix < Nx; ++ix) {
      double x_i = -xgrid/2 + spacing/2 + ix * spacing;
      for (int iz = 0; iz < Nz; ++iz) {
        double z_i = -zgrid/2 + spacing/2 + iz * spacing;

        double ex_sum = 0.0;
        double ez_sum = 0.0;

        for (int mx = 0; mx < Nx; ++mx) {
          double x_m = -xgrid/2 + spacing/2 + mx * spacing;
          for (int mz = 0; mz < Nz; ++mz) {
            double z_m = -zgrid/2 + spacing/2 + mz * spacing;

            if (ix == mx && iz == mz) continue; // skip self-contribution

            double dx = x_i - x_m;
            double dz = z_i - z_m;
            double r2 = dx*dx + dz*dz;

            double Q = chargeDensity[mx][mz] * spacing * (yMax - yMin) * spacing; // voxel charge fC

            double r3 = std::pow(r2, 1.5);
            ex_sum += k * Q * dx / r3;
            ez_sum += k * Q * dz / r3;
          }
        }
        Ex[ix][iz] = ex_sum;
        Ez[ix][iz] = ez_sum;
      }
    }  

    // Save Ex and Ez to a file in XYZ format
    std::ofstream efieldFile("electric_field.xyz");
    if (!efieldFile) {
      std::cerr << "Cannot open file for writing electric field.\n";
    } else {
      efieldFile << std::scientific << std::setprecision(6);
      const double y = 0.0; // single slice in y
      for (int ix = 0; ix < Nx; ++ix) {
        double x = -xgrid/2 + spacing/2 + ix * spacing;
        for (int iz = 0; iz < Nz; ++iz) {
            double z = -zgrid/2 + spacing/2 + iz * spacing;
            double ex = Ex[ix][iz];
            double ey = 0.0; // assume no Ey
            double ez = Ez[ix][iz];
            efieldFile << x << " " << y << " " << z << " "
                       << ex << " " << ey << " " << ez << "\n";
        }
      }
      efieldFile.close();
      std::cout << "Electric field saved to electric_field.xyz\n";
    } 

    // ------------------------- End Space Charge ------------------------------------

    // drift the positive and negative ions
    std::cout << "Negative Ions: " << drift.GetNegativeIons().size() << std::endl;
    std::cout << "Positive Ions: " << drift.GetIons().size() << std::endl;
    drift.SetTimeWindow(t, t + dt);
    drift.ResumeAvalanche(); // drift the ions

    // record recombined particles
    for (const auto& ion : drift.GetIons()) {
                if (ion.status == -9) {
                  const auto& p1 = ion.path.back();
                  ion_recombination_positions.push_back({p1.x, p1.y, p1.z});
                  recombine_num += 1;
                }
            }
    for (const auto& negion : drift.GetNegativeIons()) {
                if (negion.status == -9) {
                  const auto& p1 = negion.path.back();
                  negion_recombination_positions.push_back({p1.x, p1.y, p1.z});
                  recombine_num += 1;
                }
            }
    
    if (generateVideo && (static_cast<int>(t/dt) % videoInterval == 0)) {
      std::ostringstream ss;
      ss << std::fixed << std::setprecision(1) << (t / 1000.0);
      std::string t_str = ss.str();
        // save particle positions
        if (!aval.GetElectrons().empty()) {
            std::ofstream outfile;
            outfile.open("particle_positions/electrons/" + t_str + "us.txt", std::ios::out);
            for (const auto& electron : aval.GetElectrons()) {
                const auto& p1 = electron.path.back();
                outfile << p1.x << " " << p1.y << " " << p1.z << "\n";
            }
          outfile.close();
        }
        if (!drift.GetIons().empty()) {
            std::ofstream outfile;
            outfile.open("particle_positions/ions/" + t_str + "us.txt", std::ios::out);
            for (const auto& ion : drift.GetIons()) {
                if (ion.status == -9) {std::cout << "adding recombined ion";}
                const auto& p1 = ion.path.back();
                outfile << p1.x << " " << p1.y << " " << p1.z << "\n";
            }
          outfile.close();
        }
        
        if (!drift.GetNegativeIons().empty()) {
            std::ofstream outfile;
            outfile.open("particle_positions/negions/" + t_str + "us.txt", std::ios::out);
            for (const auto& negion : drift.GetNegativeIons()) {
                const auto& p1 = negion.path.back();
                outfile << p1.x << " " << p1.y << " " << p1.z << "\n";
            }
        outfile.close();    
        }

    }


    
    std::cout << t + dt << "ns simulated" << std::endl;
    t += dt;
    particleNum = aval.GetElectrons().size() + drift.GetIons().size() + drift.GetNegativeIons().size();
  }
  std::cout << "Recombined particles : " << recombine_num << std::endl;



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
    outfile.close();
    }
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

  std::cout << "Elapsed time: " << elapsed_ms << " ms" << std::endl;


  app.Run(kTRUE);
}