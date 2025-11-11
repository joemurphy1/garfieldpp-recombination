#include <TApplication.h>
#include <TCanvas.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <chrono>
#include <numeric>

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

// my function to integrate a Trapezoid.
double integrateTrapezoid(const double* signal, const double* time, size_t time_steps) {
    if (time_steps < 2) return 0.0; // need at least 2 points

    double integral = 0.0;

    for (size_t i = 0; i < time_steps - 1; ++i) {
        double dt = time[i + 1] - time[i];
        integral += 0.5 * (signal[i] + signal[i + 1]) * dt;
    }

    return integral;
}


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

  //Output test values
  std::cout << "Seed: " << seed << std::endl;
  std::cout << "First random number ROOT gRandom: " << gRandom->Rndm() << std::endl;
  std::cout << "First random number Garfield " << Garfield::Random::Draw() << std::endl;

  // Make a gas medium.
  

  MediumMagboltz gas;
  gas.LoadGasFile("n2_78.08_o2_20.95_ar_0.93_co2_0.04_1atm.gas");
  gas.LoadIonMobility("IonMobility_N2+_N2.txt");
  gas.LoadNegativeIonMobility("NegIonMobility_O2-_air.txt");
 
  // Make a component with analytic electric field.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  // Plate Seperation (cm)
  const double xMin = -0.15, xMax = 0.15;
  // Width in beam direction
  const double yMin = -0.15, yMax = 0.15;
  // Length? Of the plates (vertically)
  const double zMin = -0.15, zMax = 0.15;
  // Voltages
  const double vAnode = 15.;
  const double vCathode = 0.;
 // add the cathode and anode plates#]


 cmp.AddPlaneX(xMin, vCathode);
 cmp.AddPlaneX(xMax, vAnode, "anode");

   // Grid parameters.
  const int Nx = 1000;
  const int Ny = 1;
  const int Nz = 1000;
  const double spacing = 0.0001; //(cm)
  const double xgrid = Nx * spacing;
  const double zgrid = Nz * spacing;
  const double alpha = 1.72e-15; // recombination coefficient (cm^3/ns)


  ComponentGrid grid;
  grid.SetMedium(&gas);
  grid.SetMesh(Nx, Ny, Nz, -xgrid/2, xgrid/2, -xgrid/2,
                 xgrid/2, -zgrid/2, zgrid/2);
  grid.SetUniformElectricField(-vAnode/(xMax-xMin), 0., 0.); 

  // Make a sensor.
  Sensor sensor(&cmp);
  sensor.AddElectrode(&cmp, "anode");
  sensor.AddComponent(&grid);

  // option to make strip sensors for positional resolution
  const bool stripSensor = true;
  const int nStrips = 20;
  const double stripWidth = (yMax - yMin) / nStrips;
  std::vector<std::string> stripNames;
  if (stripSensor) {
    for (double y_0 = yMin; (y_0 + stripWidth) < (yMin + nStrips * stripWidth); y_0 += stripWidth ) {
      stripNames.push_back(std::to_string(y_0));
      cmp.AddStripOnPlaneX('z', xMax, y_0, y_0 + stripWidth, stripNames.back());
      sensor.AddElectrode(&cmp, stripNames.back());
    }
  }

  // Set the signal time window.
  const double tstep = 50; // monte-carlo step size (ns)
  const double tmin = -0.5 * tstep; 
  const std::size_t nbins = 50000;
  const double dt = 100.; //time between loops (dt > tstep)
  const bool stop_at_max_time = true;
  const size_t max_time = nbins*tstep;
  sensor.SetTimeWindow(tmin, tstep, nbins);
  // Set the delta reponse function.
  if (!readTransferFunction(sensor)) return 0;
  sensor.ClearSignal();

  // Set up Heed.
  TrackHeed track(&sensor);
  track.SetParticle("proton");
  track.SetEnergy(250.e6 + ProtonMass);


  // RKF integration.
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
  
  
  
  sensor.ClearSignal();
  for (std::size_t j = 0; j < nTracks; ++j) {
    double offset = randomEngine.Draw();
    double x_proton = x0 + (offset - 0.5) * 0.01; //spread over 0.1 mm
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

    // output positions and number of particles for checking
    std::cout << "Positive Ions: " << drift.GetIons().size() << std::endl;
    const auto& negions = drift.GetNegativeIons();
    if (!negions.empty()) {
      const auto& p1 = negions.front().path.back();
      std::cout << "Negative Ion at (" 
              << p1.x << ", " << p1.y << ", " << p1.z << ")\n";
    }

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

    std::cout << "Negative Ions: " << drift.GetNegativeIons().size() << std::endl;
    drift.SetTimeWindow(t, t + dt);
    drift.ResumeAvalanche(); // drift the ions

    grid.ClearFields();  // clear old densities/fields
    grid.SetUniformElectricField(-vAnode/(xMax-xMin), 0., 0.);  // maintain applied field

    const int multiplicity = 1;  // or however many charges per ion you want
    for (const auto& ion : drift.GetIons()) {
      if (!ion.path.empty()) {
        const auto& p1 = ion.path.back();
        grid.AddIon(p1.x, p1.y, p1.z, multiplicity);
      }
    }

    for (const auto& ion : drift.GetIons()) {
                if (ion.status == -9) {
                  std::cout << "Recombined Ion" << std::endl;
                  recombine_num += 1;
                }
            }
    for (const auto& negion : drift.GetNegativeIons()) {
                if (negion.status == -9) {
                  std::cout << "Recombined Negative Ion" << std::endl;
                  recombine_num += 1;
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
  sensor.IntegrateSignal("anode");
  double total_charge = sensor.GetSignal("anode", nbins - 1);
  std::cout << "Total collected charge: " << total_charge << " fC" << std::endl;
  std::cout << "Corresponding number of electrons: " << total_charge / ElementaryCharge << std::endl;
}


// option to integrate signal with trapezoidal rule doesn't save over signal
bool integrateSignalWithTrapz = false;
if (integrateSignalWithTrapz) {
  double tstart_new, tstep_new;
  size_t nsteps_new;

  sensor.GetTimeWindow(tstart_new, tstep_new, nsteps_new);
  std::vector<double> signal(nsteps_new);
  std::vector<double> time(nsteps_new);

  for (size_t step_number = 0; step_number < nsteps_new; ++step_number) {
      signal[step_number] = sensor.GetSignal("anode", step_number);
      time[step_number] = tstart_new + (step_number * tstep_new);
  }


  double integral = integrateTrapezoid(signal.data(), time.data(), nsteps_new);
  std::cout << "Integrates to: " << integral << std::endl;
}


  bool saveSignal = true;
  if (saveSignal) { std::ofstream outfile;
    outfile.open("signal.txt", std::ios::out);
    for (unsigned int i = 0; i < nbins; ++i) {
    const double t = (i + 0.5) * tstep;
    const double f = sensor.GetSignal("anode", i);
    const double fe = sensor.GetElectronSignal("anode", i);
    const double fh = sensor.GetIonSignal("anode", i);
    outfile << t << " " << f << " " << fe << " " << fh << "\n";
    outfile.close();
    }
    if (stripSensor) {
      std::cout << "Saving strip signals to file." << std::endl;
      sensor.IntegrateSignals();
      outfile.open("strip_signals.txt", std::ios::out);
      outfile << "Strip_y_position(cm) total_charge(fC)\n";
      for (const auto& name : stripNames) {
        const double strip_total_charge = sensor.GetSignal(name, nbins - 1);
        outfile << name << " " << strip_total_charge <<"\n";
      }
    }
  }


  
  if (plotSignal) { 
    sensor.PlotSignal("anode", cS);
  }

  // timer
  auto time_end = std::chrono::steady_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();

  std::cout << "Elapsed time: " << elapsed_ms << " ms" << std::endl;


  app.Run(kTRUE);
}