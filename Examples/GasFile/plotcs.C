#include <iostream>

#include "TCanvas.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  const double pressure = AtmosphericPressure;
  const double temperature = 293.15;

  // Setup the gas component 1.
  MediumMagboltz gascmp1("N2", 100.);
  gascmp1.SetTemperature(temperature);
  gascmp1.SetPressure(pressure);

  TCanvas * c1 = new TCanvas("c1", "c1", 800,800);
  gascmp1.SetCanvas(c1);
  gascmp1.PlotElectronCrossSections();
  c1->SaveAs("DryAir_CrossSections_N2.png");

  // Setup the gas component 2.
  MediumMagboltz gascmp2("O2", 100.);
  gascmp2.SetTemperature(temperature);
  gascmp2.SetPressure(pressure);

  TCanvas * c2 = new TCanvas("c2", "c2", 800,800);
  gascmp2.SetCanvas(c2);
  gascmp2.PlotElectronCrossSections();
  c2->SaveAs("DryAir_CrossSections_O2.png");
  

  
  // Setup the gas mixture.
  MediumMagboltz gasmix("N2", 80.1, "O2", 19.9);
  gasmix.SetTemperature(temperature);
  gasmix.SetPressure(pressure);
  gasmix.Initialise();
  
  TCanvas * c3 = new TCanvas("c3", "c3", 800,800);
  gasmix.SetCanvas(c3);
  gasmix.PlotElectronCollisionRates();
  c3->SaveAs("DryAir_CollisionRates.png");
  
}
