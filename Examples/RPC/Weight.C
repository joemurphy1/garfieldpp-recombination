#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>

#include "Garfield/AvalancheGrid.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewSignal.hh"

#define LOG(x) std::cout << x << std::endl

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  const bool debug = true;
  constexpr bool plotSignal = true;
  constexpr bool plotDrift = false;
    
  // The geometry of the RPC.
  const double gap = 0.03;  // [cm]
  const double thickness = 0.2;  // [cm]
  const double eps = 8.;  // [1]
  const double voltage = -3e3;  // [V]
  const double sigma = 1e-12; // [S/m]

  ComponentParallelPlate* RPC = new ComponentParallelPlate();
  RPC->Setup(gap, thickness, eps, voltage, sigma);

  // Adding a readout structure.
  const std::string label = "ReadoutPlane";
  RPC->AddPlane(label);

  // Setup the gas, but one can also use a gasfile.
  MediumMagboltz gas;
  gas.LoadGasFile("TimingRPCGas.gas"); // c2h2f4/ic4h10/sf6 85/5/10.
  gas.Initialise(true);

  // Setting the drift medium.
  RPC->SetMedium(&gas);

  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(RPC);
  sensor.AddElectrode(RPC, label);

  // The resistive layer will make the weighting potential time-dependent.
  sensor.EnableDelayedSignal();

  // Set the time bins.
  const unsigned int nTimeBins = 200;
  const double tmin = 0.;
  const double tmax = 4;
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);

  // Set the times the delayed signals will be calculated.
  std::vector<double> times;
  for (int i = 0; i < nTimeBins; i++) {
    times.push_back(tmin + tstep / 2 + i * tstep);
  }
  //sensor.SetDelayedSignalTimes(times);

  // Create the AvalancheMicroscopic.
  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);
  aval.UseWeightingPotential();
  aval.EnableSignalCalculation();

  // Set time window where the calculations will be done microscopically.
  const double tMaxWindow = 2;
  aval.SetTimeWindow(0., tMaxWindow);

  // Create the AvalancheGrid for grid-based avalanche calculations that are
  // suited for constant drift fields. This class will take over the
  // calculations of the microscopic class after the set time-window.
  AvalancheGrid avalgrid;
  avalgrid.SetSensor(&sensor);
  avalgrid.SetAvalancheMicroscopic(&aval);
    
    
  avalgrid.SetGrid(-0.05, 0.05, 5,-0.05, 0.05, 5,0, gap, 1000);
  //avalgrid.EnableDebugging();

  // Preparing the plotting of the induced charge and signal of the electrode
  // readout.
  ViewSignal* signalView = nullptr;
  TCanvas* cSignal = nullptr;
  if (plotSignal) {
    cSignal = new TCanvas("cSignal", "", 600, 600);
    signalView = new ViewSignal();
    signalView->SetCanvas(cSignal);
    signalView->SetSensor(&sensor);
  }

  ViewSignal* chargeView = nullptr;
  TCanvas* cCharge = nullptr;

  if (plotSignal) {
    cCharge = new TCanvas("cCharge", "", 600, 600);
    chargeView = new ViewSignal();
    chargeView->SetCanvas(cCharge);
    chargeView->SetSensor(&sensor);
  }

  // Set up Heed.
  TrackHeed track;
  track.SetSensor(&sensor);
  // Set the particle type and momentum [eV/c].
    track.SetParticle("pion");
    track.SetMomentum(7.e9);

  ViewDrift* driftView = nullptr;
  TCanvas* cDrift = nullptr;
  if (plotDrift) {
    cDrift = new TCanvas("cDrift", "", 600, 600);
    driftView = new ViewDrift();
    driftView->SetArea(-0.2, -0.2, 0, 0.2, 0.2, gap);
    driftView->SetPlane(0, -1, 0, 0, 0, 0);
    driftView->SetCanvas(cDrift);
    track.EnablePlotting(driftView);
  }

  // Setting the timer for the running time of the algorithm.
  std::clock_t start = std::clock();

  // Simulate a charged-particle track.
  track.NewTrack(0, 0, gap, 0, 0, 0, -1);
  // Retrieve the clusters along the track.
  for (const auto& cluster : track.GetClusters()) {
    // Loop over the electrons in the cluster.
    for (const auto& electron : cluster.electrons) {
      // Simulate the electron and hole drift lines.
      if (plotDrift) aval.EnablePlotting(driftView);
      aval.AvalancheElectron(electron.x, electron.y, electron.z, 
                             electron.t, 0.1, 0., 0., 0.);
      // Stops calculation after tMaxWindow ns.
      avalgrid.GetElectronsFromAvalancheMicroscopic();
    }
  }

  // Start grid based avalanche calculations starting from where the microsocpic
  // calculations stoped.
  avalgrid.StartGridAvalanche();
  // Stop timer.
  double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;

  LOG("Script: "
      << "Electrons have drifted. It took " << duration << "s to run.");
    
  if (plotDrift) {
    constexpr bool twod = true;
    driftView->Plot(twod);
    cDrift->Update();
    gSystem->ProcessEvents();
  }
    
  if (plotSignal) {
    // Plot signals
    signalView->Plot(label, true);
    cSignal->Update();
    gSystem->ProcessEvents();
      
    sensor.ExportSignal(label,"SignalNaN");
    // Plot induced charge
    chargeView->Plot(label, false);
    cCharge->Update();
    gSystem->ProcessEvents();
    // Export induced current data as an csv file.
    sensor.ExportSignal(label,"ChargeNaN");
  }
    
  app.Run(kTRUE);
}
