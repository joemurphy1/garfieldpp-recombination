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
    std::vector<double> v_ion;
    std::vector<double> v_negion;

    for (double E_x = 0; E_x < 50000; E_x += 100) {
        double vx_ion, vy_ion, vz_ion;
        double vx_negion, vy_negion, vz_negion;
        // Query the medium directly for ion velocities at the given E-field
        gas.IonVelocity(E_x, 0, 0, 0, 0, 0, vx_ion, vy_ion, vz_ion); // cm/ns
        gas.NegativeIonVelocity(E_x, 0, 0, 0, 0, 0, vx_negion, vy_negion, vz_negion); // cm/ns

        Ex.push_back(E_x);
        v_ion.push_back(vx_ion);
        v_negion.push_back(vx_negion);
    }

    std::ofstream outfile("E_field_Ion_velocity.txt");
    outfile << "Electric Field [V/cm] | N2+ Velocity [cm/ns] | O2- Velocity [cm/ns]\n";
    for (size_t i = 0; i < Ex.size(); ++i) {
        outfile << Ex[i] << " " << v_ion[i] << " " << v_negion[i] << "\n";
    }
    outfile.close();

    // If you want to show a GUI or keep the application running, uncomment:
    // app.Run();

    return 0;
}




