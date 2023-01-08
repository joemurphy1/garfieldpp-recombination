#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackBichsel.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  randomEngine.Seed(123456);
  TApplication app("app", &argc, argv);
  SetDefaultStyle();

  // Histograms
  TH1::StatOverflows(true); 
  TH1F hEdep("hEdep", "Energy Loss", 100, 0., 10.);

  MediumSilicon si;

  constexpr double width = 10.e-4;

  // Make a component
  ComponentConstant cmp;
  cmp.SetArea(0., -10., -10., width, 10., 10.);
  cmp.SetMedium(&si);
  cmp.SetElectricField(100., 0., 0.);

  // Make a sensor
  Sensor sensor;
  sensor.AddComponent(&cmp);

  TrackBichsel track;
  track.EnableDebugging();
  track.SetSensor(&sensor);
  track.SetParticle("pi");
  track.SetBetaGamma(10.);
  track.Initialise();
  track.ComputeCrossSection();
  track.DisableDebugging();
  const int nEvents = 10000;
  for (int i = 0; i < nEvents; ++i) {
    if (i % 1000 == 0) std::cout << i << "/" << nEvents << "\n";
    // Initial position and direction 
    double x0 = 0., y0 = 0., z0 = 0., t0 = 0.;
    double dx0 = 1., dy0 = 0., dz0 = 0.; 
    track.NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
    // Cluster coordinates
    double xc = 0., yc = 0., zc = 0., tc = 0.;
    // Number of electrons produced in a collision
    int nc = 0;
    // Energy loss in a collision
    double ec = 0.;
    // Dummy variable (not used at present)
    double extra = 0.;
    // Total energy loss along the track
    double esum = 0.;
    // Total number of electrons produced along the track
    int nsum = 0;
    // Loop over the clusters.
    while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      esum += ec;
    }
    hEdep.Fill(esum * 1.e-3);
  }
 
  TCanvas c1;
  hEdep.GetXaxis()->SetTitle("energy loss [keV]");
  hEdep.Draw();
  c1.SaveAs("edep.pdf");

  app.Run(true); 

}
