//
// Created by Dario Stocco (stoccod@ethz.ch) on 02.08.2023.
//
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <numeric>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TSystem.h>

#include "Garfield/AvalancheGridSpaceCharge.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewSignal.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {
  std::cout << "Start RPC Space Charge Example\n";

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  double voltage = 9000;

  MediumMagboltz gas;
  gas.EnableAutoEnergyLimit(false);
  gas.SetMaxElectronEnergy(100.);
  std::string path = "CO2_C2H2F4_isoC4H10_SF6_30_64.5_4.5_1_T_293_P_723.8.gas";
  gas.LoadGasFile(path);
  gas.Initialise(true);

  // Dimensions of RPC
  double d_bakelite = 0.2; // (cm)
  double d_pet = 0.02;
  double d_gas = 0.2;
  std::vector<double> layers = {d_pet,
                                d_bakelite,
                                d_gas,
                                d_bakelite,
                                d_pet};
  double y_mid = d_pet + d_bakelite + d_gas / 2;

  double e_bakelite = 8.;
  double e_pet = 3.5;
  double e_gas = 1.;
  std::vector<double> eps = {e_pet,
                             e_bakelite,
                             e_gas,
                             e_bakelite,
                             e_pet};
  // ComponentParallelPlate
  ComponentParallelPlate cmp;
  cmp.Setup(int(layers.size()), eps, layers, voltage, {});
  std::string label = "readout";
  cmp.AddPlane(label);
  cmp.SetMedium(&gas);

  // Sensor
  Sensor sens;
  sens.AddComponent(&cmp);
  sens.AddElectrode(&cmp, label);
  sens.SetTimeWindow(0, (25. - 0) / 200., 200);

  // AvalancheGridSpace Charge
  AvalancheGridSpaceCharge avalsc;
  avalsc.EnableDebugging();
  avalsc.EnableDiffusion(true);
  avalsc.EnableStickyAnode(true);
  avalsc.EnableAdaptiveTimeStepping(true);
  avalsc.SetStopAtK(true);
  // disable space charge calculation
  avalsc.EnableSpaceChargeEffect(false);
  // set sensor and grid
  avalsc.SetSensor(&sens);
  avalsc.Set2dGrid(y_mid - 0.5 * d_gas + 1.e-8, 
                   y_mid + 0.5 * d_gas - 1.e-8, 400, 0.05, 100);

  // Avalanche electron
  // place 1000 electrons in the middle of the gas gap
  avalsc.AvalancheElectron(0., y_mid, 0., 0., 1000.);
  avalsc.StartGridAvalanche();
  // export grid
  avalsc.ExportGrid("my_rpc_grid");

  // view recorded signals from plane electrode
  ViewSignal *signal_view = new ViewSignal(&sens);
  TCanvas *c_signal = new TCanvas(label.c_str(), label.c_str(), 600, 600);
  signal_view->SetCanvas(c_signal);
  signal_view->PlotSignal(label);
  c_signal->SetTitle(label.c_str());
  gSystem->ProcessEvents();

  app.Run(true);
  return 0;
}
