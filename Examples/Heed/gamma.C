#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "Garfield/TrackHeed.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  // Make a gas medium.
  MediumMagboltz gas("Ar", 70., "CO2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(AtmosphericPressure);
  
  constexpr double r = 100.;
  // SolidTube tube(0, 0, 0, r, 10000.);
  SolidTube tube(0, 0, 0, r, 10.);

  // Combine gas and box to a simple geometry.
  GeometrySimple geo;
  geo.AddSolid(&tube, &gas);

  // Make a component with constant electric field.
  ComponentConstant field;
  field.SetGeometry(&geo);
  field.SetElectricField(0., 0., 500.); 

  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&field);
  
  // Use Heed for simulating the photon absorption.
  TrackHeed track;
  track.SetSensor(&sensor);
  track.DisableDeltaElectronTransport();

  TH1F hPhi("hPhi", "phi", 360, -Pi, Pi);
  TH1F hTheta("hTheta", "ctheta", 100, -1, 1.);

  // const unsigned int nEvents = 500000;
  const unsigned int nEvents = 1000;
  for (unsigned int i = 0; i < nEvents; ++i) {
    if (i % 1000 == 0) std::cout << i << "/" << nEvents << "\n";
    // Initial coordinates of the photon.
    const double egamma = 100.e3;
    int ne = 0;
    // track.TransportPhoton(0., 0., 0., 0., egamma, 0., 0., 1., ne);
    track.NewTrack(0., 0., 0., 0., 0., 0., 1.);
    double xc, yc, zc, tc, ec, extra;
    while (track.GetCluster(xc, yc, zc, tc, ne, ec, extra)) {
      for (int j = 0; j < ne; ++j) {
        double xe, ye, ze, te, ee, dxe, dye, dze;
        track.GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);
        if (j == 0) {
          // std::cout << "  (" << dxe << ", " << dye << ", " << dze << ")\n";
        // } else {
          hTheta.Fill(dze);
          const double theta = acos(dze);
          const double phi = acos(dxe / sin(theta));
          hPhi.Fill(phi); 
        }
      }
    }
  }

  TCanvas c("c", "", 600, 600);
  c.Divide(2, 1);
  c.cd(1);
  hPhi.Draw();
  c.cd(2);
  hTheta.Draw();
  c.Update();
  app.Run(true);
}
