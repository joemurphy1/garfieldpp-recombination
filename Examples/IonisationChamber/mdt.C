#include <TApplication.h>
#include <TCanvas.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <chrono>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/RandomEngine.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Random.hh"
#include "Garfield/FundamentalConstants.hh"


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
  gas.LoadGasFile("n2_79_o2_21_1atm.gas");
  gas.LoadIonMobility("IonMobility_N2+_N2.txt");
  gas.LoadNegativeIonMobility("NegIonMobility_O2-_air.txt");


  // Make a component with analytic electric field.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  // Plate Seperation (cm I think)
  const double xMin = -15, xMax = 15;
  // Width in beam direction
  const double yMin = -4, yMax = 4;
  // Length? Of the plates (vertically)
  const double zMin = -15, zMax = 15;
  // Voltages
  const double vAnode = 1500.;
  const double vCathode = 0.;
 // add the cathode and anode plates#]


 cmp.AddPlaneX(xMin, vCathode);
 cmp.AddPlaneX(xMax, vAnode, "anode");

  // Make a sensor.
  Sensor sensor(&cmp);
  sensor.AddElectrode(&cmp, "anode");
  // Set the signal time window.
  const double tstep = 0.5;
  const double tmin = -0.5 * tstep;
  const std::size_t nbins = 10000;
  sensor.SetTimeWindow(tmin, tstep, nbins);
  // Set the delta reponse function.
  if (!readTransferFunction(sensor)) return 0;
  sensor.ClearSignal();

  // Set up Heed.
  TrackHeed track(&sensor);
  track.SetParticle("proton");
  track.SetEnergy(250.e6 + ProtonMass);


  // RKF integration.
  DriftLineRKF drift(&sensor);
  drift.SetGainFluctuationsPolya(0., 20000.);
  drift.EnableNegativeIonTail(true);

  TCanvas* cD = nullptr;
  ViewDrift driftView;
  constexpr bool plotDrift = true;
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
  const std::size_t nTracks = 1;

  
  sensor.ClearSignal();
  for (std::size_t j = 0; j < nTracks; ++j) {
    track.NewTrack(x0, y0, z0, 0, 0, 1, 0);
    for (const auto& cluster : track.GetClusters()) {
      //remove clusters that are unphysically out of the detector
      if (cluster.y < yMin || cluster.y > yMax) {
        continue;
      }

      for (const auto& electron : cluster.electrons) {
        drift.DriftElectron(electron.x, electron.y, electron.z, electron.t);
      }
      for (const auto& negativeIon : cluster.ions) {
        drift.DriftNegativeIon(negativeIon.x, negativeIon.y, negativeIon.z, negativeIon.t);
      }
    }
    if (plotDrift) {
      cD->Clear();
      cmp.PlotCell(cD);
      constexpr bool twod = true;
      constexpr bool drawaxis = false;
      driftView.Plot(twod, drawaxis);
    }
    
  }

  //sensor.ConvoluteSignals();
  int nt = 0;

    // get all the data by looping through the bins.

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
  // option to display integrated signal
  bool display_integrated_signal = false;
  if (display_integrated_signal) sensor.IntegrateSignals();

  bool saveSignal = true;
  if (saveSignal) { std::ofstream outfile;
    outfile.open("signal.txt", std::ios::out);
    for (unsigned int i = 0; i < nsteps_new; ++i) {
    const double t = (i + 0.5) * tstep_new;
    const double f = sensor.GetSignal("anode", i);
    const double fe = sensor.GetElectronSignal("anode", i);
    const double fh = sensor.GetIonSignal("anode", i);
    outfile << t << " " << f << " " << fe << " " << fh << "\n";
    }
  }
  //std::cout << sensor.GetSignal("s", nsteps_new-1) << std::endl;
  //std::cout << sensor.GetSignal("tube", nsteps_new-1) << std::endl;
  if (sensor.ComputeThresholdCrossings(-2., "anode", nt)) {
    if (plotSignal) { 
      sensor.PlotSignal("anode", cS);
    }
  }

  // timer
  auto time_end = std::chrono::steady_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();

  std::cout << "Elapsed time: " << elapsed_ms << " ms" << std::endl;


  app.Run(kTRUE);
}