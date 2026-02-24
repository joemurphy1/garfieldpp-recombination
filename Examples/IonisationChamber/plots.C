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

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);

    MediumMagboltz gas;
    gas.LoadGasFile("n2_78.08_o2_20.95_ar_0.93_co2_0.04_1atm.gas");
    gas.LoadIonMobility("IonMobility_N2+_N2.txt");
    gas.LoadNegativeIonMobility("NegIonMobility_O2-_air.txt");

    ComponentAnalyticField cmp;
    cmp.SetMedium(&gas);

    std::vector<double> Ex;
    std::vector<double> v;
    std::vector<double> t;

    for (double E_x = 0; E_x < 50000; E_x += 1) {
        double vx, vy, vz;
        double eta;
        // Query the medium directly for ion velocities at the given E-field
        gas.ElectronVelocity(E_x, 0, 0, 0, 0, 0, vx, vy, vz); // cm/ns
        gas.ElectronAttachment(E_x, 0, 0, 0, 0, 0, eta); // cm/ns

        Ex.push_back(E_x);
        v.push_back(vx);
        t.push_back(1/(eta*vx));
    }

    std::ofstream outfile("E_field_Ion_velocity.txt");
    outfile << "Electric Field [V/cm] | Electron Velocity [cm/ns] | Attachment lifetime [ns]\n";
    for (size_t i = 0; i < Ex.size(); ++i) {
        outfile << Ex[i] << " " << v[i] << " " << t[i] << "\n";
    }
    outfile.close();

    // If you want to show a GUI or keep the application running, uncomment:
    // app.Run();

    return 0;
}




