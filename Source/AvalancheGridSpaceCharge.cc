#include "Garfield/AvalancheGridSpaceCharge.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/Exceptions.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Medium.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"

namespace {

double Mag(const double x, const double y) { return std::sqrt(x * x + y * y); }

// Get size of avalanche when going from x to x+dx in Monte Carlo fashion
void GetAvalancheSizeFromStep(double dx, const long nEIn,
                              const double alpha, const double eta,
                              long &nEOut, double &nPOut, double &nNOut) {
  // Monte Carlo Avalanche gain per travelled distance dx (cm)
  nEOut = 0;
  nPOut = 0;
  nNOut = 0;

  if (std::abs(alpha - eta) < 1.e-8 && alpha > 1.e-8) {
    // alpha == eta
    if (nEIn < 1000L) {
      // Condition to which the random number will be compared.
      // If the number is smaller than the condition, nothing happens.
      // Otherwise, the single electron will be attached or retrieve
      // additional electrons from the gas.
      const double p = alpha * dx / (1. + alpha * dx);
      const double q = 1. + alpha * dx;
      const double f = 1. / log(p);
      // Running over all electrons in the avalanche.
      for (long i = 0; i < nEIn; i++) {
        // Draw a random number from the uniform distribution (0,1).
        const double s = Garfield::RndmUniformPos();
        // We (wrongly) assume if s >= prob only pos ions are created
        // else 1 neg ion.
        if (s >= p) {
          // Deviation/improvement wrt Lippmann?
          nEOut += (long)(log((1. - s) * q) * f);
        } else {
          nNOut += 1;
        }
      }
      // charge conservation
      nPOut = (nEOut - nEIn) + nNOut;

    } else {
      // Central limit theorem.
      const double sigma = sqrt(2 * alpha * dx * nEIn);
      nEOut = (long)Garfield::RndmGaussian(nEIn, sigma);

      // boundary conditions (alpha dx ElectronIn = dPosOut),
      //  the procedure guarantees positive values
      //  and net Nion = nPos - nNeg = nOut - nIn
      if (nEOut <= 0) nEOut = 0;  //< unphysical
      if (nEOut >= nEIn) {
        nNOut = std::expm1(eta * dx) * nEIn;         //< >= 0
        nPOut = (nEOut - nEIn) + nNOut;  //< >= 0
      } else {
        nPOut = std::expm1(alpha * dx) * nEIn;       //< >= 0
        nNOut = nPOut - (nEOut - nEIn);  //< >= 0
      }
    }
  } else if (alpha < 1.e-8 && eta > 1.e-8) {
    // alpha == 0, only attachment possible
    if (nEIn < 1000L) {
      const double p = exp(-eta * dx);
      for (long i = 0; i < nEIn; i++) {
        // Draw a random number from the uniform distribution (0,1).
        const double s = Garfield::RndmUniformPos();
        if (s >= p) {
          nNOut += 1;
        } else {
          nEOut += 1;
        }
      }
    } else {
      // Central limit theorem.
      const double ndx = exp(-eta * dx);
      const double mu = nEIn * ndx;
      const double sigma = std::sqrt(mu * ndx * (ndx - 1.));
      nEOut = (long)Garfield::RndmGaussian(mu, sigma);

      // boundary conditions
      if (nEOut <= 0) {
        // Unphysical
        nEOut = 0;
      } else if (nEOut > nEIn) {
        // Unphysical with alpha = 0
        nEOut = nEIn;
      }
      // Charge conservation
      nNOut = -(nEOut - nEIn);
    }
  } else {
    // alpha != 0 =! eta
    const double k = eta / alpha;
    const double ndx = exp((alpha - eta) * dx);

    if (nEIn < 1000L) {
      // Condition to which the random number will be compared.
      // If the number is smaller than the condition, nothing happens.
      // Otherwise, the single electron will be attached or retrieve
      // additional electrons from the gas.
      const double p = k * (ndx - 1.) / (ndx - k);
      const double q = (ndx - k) / (ndx * (1. - k));
      const double f = 1. / log(1. - (1. - k) / (ndx - k));
      // Running over all electrons in the avalanche.
      for (long i = 0; i < nEIn; i++) {
        // Draw a random number from the uniform distribution (0,1).
        const double s = Garfield::RndmUniformPos();
        if (s >= p) {
          // deviation/improvement wrt Lippmann?
          nEOut += (long)(1. + log((1. - s) * q) * f);
        } else {
          nNOut += 1;
        }
      }
      // charge conservation
      nPOut = (nEOut - nEIn) + nNOut;

    } else {
      // Central limit theorem.
      const double mu = nEIn * ndx;
      const double sigma = std::sqrt(mu * (1. + k) * (ndx - 1.) / (1. - k));
      nEOut = (long)Garfield::RndmGaussian(mu, sigma);

      // boundary conditions (alpha dx ElectronIn = dPosOut),
      //  the procedure guarantees positive values
      //  and netto Nion = nPos - nNeg = nOut - nIn
      if (nEOut <= 0) nEOut = 0;  //< unphysical

      // either the above has not been executed or nPOut was not positive
      if (nEOut >= nEIn) {
        nNOut = std::expm1(eta * dx) * nEIn;
        // nNOut = eta / (alpha - eta) * (nEOut - nEIn);
        nPOut = (nEOut - nEIn) + nNOut;
      } else {
        nPOut = std::expm1(alpha * dx) * nEIn;
        // nPOut = alpha / (alpha - eta) * (nEOut - nEIn);
        nNOut = nPOut - (nEOut - nEIn);
      }
    }
  }
}

// Get mean size of avalanche when going from x to x+dx
void GetMeanAvalancheSizeFromStep(double dx, const long nEIn,
                                  const double alpha, const double eta,
                                  long &nEOut, double &nPOut,
                                  double &nNOut) {
  // Mean size gain
  nPOut = 0;
  nNOut = 0;
  const double ndx = exp((alpha - eta) * dx);
  // for electrons
  nEOut = nEIn * ndx;
  if (nEOut <= 0) nEOut = 0;  //< unphysical

  // either the above has not been executed or nPOut was not positive
  if (nEOut >= nEIn) {
    nNOut = std::expm1(eta * dx) * nEIn;         //< >= 0
    nPOut = (nEOut - nEIn) + nNOut;  //< >= 0
  } else {
    nPOut = std::expm1(alpha * dx) * nEIn;       //< >= 0
    nNOut = nPOut - (nEOut - nEIn);  //< >= 0
  }
}

}  // namespace

namespace Garfield {

AvalancheGridSpaceCharge::AvalancheGridSpaceCharge() {
  m_zGrid.reserve(5000);
  m_rGrid.reserve(1000);
}

AvalancheGridSpaceCharge::AvalancheGridSpaceCharge(Sensor *sensor)
    : AvalancheGridSpaceCharge() {
  SetSensor(sensor);
}

void AvalancheGridSpaceCharge::SetSensor(Sensor *sensor) {
  if (!sensor) throw Exception("Sensor can't be nullptr");
  m_sensor = sensor;
  // Determine if one of the components is a parallel-plate one.
  m_pp = nullptr;
  const size_t nCmp = m_sensor->GetNumberOfComponents();
  for (size_t i = 0; i < nCmp; i++) {
    m_pp = dynamic_cast<ComponentParallelPlate *>(m_sensor->GetComponent(i));
    if (m_pp) break;
  }
}

int AvalancheGridSpaceCharge::GetGasGapNumber(int layerIndex) {
  auto it =
      std::find(m_vIndexGasGaps.begin(), m_vIndexGasGaps.end(), layerIndex);
  return (it != m_vIndexGasGaps.end())
             ? std::distance(m_vIndexGasGaps.begin(), it)
             : -1;
}

void AvalancheGridSpaceCharge::Reset() {
  m_time = 0.;
  m_time0 = 0.;
  m_dt = 0.;
  m_nTotElectron = 0;
  m_nTotPosIons = 0;

  m_vCoNGasLayer.resize(0);
  m_vElectrons.resize(0);
  m_vNElectronEvolution.resize(0);
  m_grid.clear();
  m_vYPointInGasGap.resize(0);
  m_vIndexGasGaps = {0};
  m_ezBkg = {0};
  m_vSaturatedGaps.resize(0);

  m_bDriftAvalanche = false;
  m_bImportAvalanche = false;
  m_bPreparedImportAvalanche = false;
  m_bFieldK = false;

  std::cout << m_className << "::Reset: Instance reset, ready to use again.\n";
}

void AvalancheGridSpaceCharge::Set2dGrid(const double zmin, const double zmax,
                                         const int zsteps, const double rmax,
                                         const int rsteps) {
  if (zmin >= zmax || zsteps <= 0 || 0 >= rmax || rsteps <= 0) {
    throw Exception("Grid is not properly defined");
  }

  // set z grid
  m_zSteps = zsteps;
  m_zStepSize = (zmax - zmin) / zsteps;
  m_zInvStep = 1. / m_zStepSize;
  for (int i = 0; i < zsteps + 1; i++) {
    m_zGrid.push_back(zmin + i * m_zStepSize);
  }

  // set r grid
  m_rSteps = rsteps;
  m_rStepSize = rmax / rsteps;
  m_rInvStep = 1. / m_rStepSize;
  for (int i = 0; i < rsteps + 1; i++) {  
    m_rGrid.push_back(0 + i * m_rStepSize);
  }

  if (m_bDebug) {
    std::cout << m_className << "::Set2dGrid: Grid created:\n"
              << "       z range = (" << zmin << "," << zmax << ").\n"
              << "       r range = (" << 0 << "," << rmax << ").\n";
  }
}

void AvalancheGridSpaceCharge::SetFieldCalculation(const std::string& option,
                                                   const int nof_approx) {
  std::string opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), toupper); 
  if (opt == "COULOMB") {
    m_fieldOption = FieldOption::Coulomb;
  } else if (opt == "MIRROR") {
    m_fieldOption = FieldOption::Mirror;
  } else {
    std::cerr << m_className << "::SetFieldCalculation: Unknown option "
              << option << ".\n";
  }
  m_iFieldApprox = nof_approx;
}

void AvalancheGridSpaceCharge::AddElectrons(AvalancheMicroscopic *avmc) {
  if (!avmc) return;

  if (!m_bImportAvalanche) {
    m_bImportAvalanche = true;

    // resize the electrons according to # gap
    if (m_pp) m_pp->IndexOfGasGaps(m_vIndexGasGaps);
    m_vElectrons.resize(m_vIndexGasGaps.size());
  }

  if (m_bDebug) std::cout << m_className << "::AddElectrons:\n";
  for (const auto &electron : avmc->GetElectrons()) {
    // Skip electrons that don't have status code "outside time window".
    if (electron.status != StatusOutsideTimeWindow) {
      if (m_bDebug) { 
        std::cout << "    Skipping electron with status " 
                  << electron.status << ".\n";
      }
      continue;
    }
    int k = 0;
    if (m_pp) {
      const int ind = m_pp->GetLayer(electron.path.back().y);
      if (ind < 0) {
        std::cerr << m_className
                  << "::AddElectrons: Electron outside component.\n";
        continue;
      }
      k = GetGasGapNumber(ind);
      if (k == -1) {
        std::cerr << m_className << "AddElectrons:\n"
                  << "    Electron is not in a gas gap, continue.\n";
        continue;
      }
    }
    // Add the electron to a vector.
    Point pt{};
    pt.x = electron.path.back().x;
    pt.y = electron.path.back().y;
    pt.z = electron.path.back().z;
    pt.t = electron.path.back().t;
    m_vElectrons[k].push_back(std::move(pt));

    if (m_bDebug)
      std::cout << m_className << "::AddElectrons: Electron added, y: "
                << electron.path.back().y << " and gas gap: " << k + 1 << "\n";
  }
}

void AvalancheGridSpaceCharge::AddElectron(const double x, const double y,
                                           const double z, const double t,
                                           const int n) {
  int gasGap = 0;
  // check if avalanche electron in a gas gap
  if (m_pp) {
    int ind;
    double eps = -1;
    if (!m_pp->GetLayer(y, ind, eps) && eps != 1.) {
      std::cerr << m_className << "AddElectron: Electron is not in a gas gap.";
      return;
    }
    // determine indices of gas gaps
    m_pp->IndexOfGasGaps(m_vIndexGasGaps);
    gasGap = GetGasGapNumber(ind);
  } else {
    // put the y-electron coord as reference
    m_vYPointInGasGap.push_back(y);
  }

  if (!m_bDriftAvalanche) {
    m_bDriftAvalanche = true;
  }

  if (m_time == 0 && m_time != t && m_bDebug)
    std::cerr << m_className
              << "::AddElectron: Overwriting start time of avalanche for t "
                 "= 0 to "
              << t << ".\n";

  m_time = t;
  m_time0 = t;

  // prepare the CoN
  for (int i = 0; i < (int)m_vIndexGasGaps.size(); i++) {
    if (i != gasGap) {
      // HS: not sure it's a good idea to initialize with NAN...
      m_vCoNGasLayer.push_back({NAN, NAN, NAN});
    } else if (i == gasGap) {
      m_vCoNGasLayer.push_back({x, y, z});
    }
  }

  if (m_vCoNGasLayer.size() == 0) {
    std::cerr << m_className << "::AddElectron: Could not determine center.\n";
  }
  Prepare2dMesh();
  if (SnapTo2dGrid(x, y, z, n, gasGap) && m_bDebug)
    std::cerr << m_className
              << "::AddElectron: Electron added at (t, x, y, z) =  (" << t
              << ", " << x << ", " << y << ", " << z << ").\n";
}

void AvalancheGridSpaceCharge::AddExtraElectron(double y, int n) {
  if (!m_bDriftAvalanche) {
    std::cerr << m_className << "::AddExtraElectron: First use AddElectron.\n";
    return;
  }

  // check if avalanche electron in a gas gap
  int gasGap = 0;
  if (m_pp) {
    int ind;
    double eps = -1;
    if (!m_pp->GetLayer(y, ind, eps) && eps != 1.) {
      std::cerr << m_className
                << "AddExtraElectron: Electron is not in a gas gap.";
      return;
    }
    gasGap = GetGasGapNumber(ind);
  }

  // nothing yet in this gas gap
  if (std::isnan(m_vCoNGasLayer[gasGap][0])) {
    m_vCoNGasLayer[gasGap] = {0, y, 0};
  }
  if (SnapTo2dGrid(m_vCoNGasLayer[gasGap][0], y, m_vCoNGasLayer[gasGap][2], n,
                   gasGap) &&
      m_bDebug)
    std::cout << m_className << "::AddExtraElectron: "
              << "Electron added at (t, x, y, z) =  (" << m_time << ", "
              << m_vCoNGasLayer[gasGap][0] << ", " << y << ", "
              << m_vCoNGasLayer[gasGap][2] << ").\n";
}

void AvalancheGridSpaceCharge::StartGridAvalanche(double dtime) {
  // avalanche the electrons until a certain delta-time OR there are no
  // electrons left in the gap
  if ((!m_bImportAvalanche && !m_bDriftAvalanche) || !m_sensor) return;

  // prepare the imported avalanche
  if (m_bImportAvalanche && !m_bPreparedImportAvalanche) {
    PrepareElectronsFromMicroscopicAvalanche();
    if (m_bDebug)
      std::cerr << m_className
                << "::StartGridAvalanche: Microscopic electrons successfully "
                   "prepared.\n";
    m_bPreparedImportAvalanche = true;
  }

  // check if electrons are on grid
  if (m_nTotElectron <= 0) {
    std::cerr << m_className << "::StartGridAvalanche: Cancelled "
              << m_nTotElectron << " electrons on grid.\n";
    return;
  }

  if (dtime == -1) {
    // transport until no electrons in gap or an error is reached
    while (true) {
      // Transport the nodes (returns false if 0 electrons in gap)
      if (!TransportTimeStep()) {
        break;
      }
    }
  } else {
    double tStart = m_time;

    while (m_time + m_dt - tStart < dtime) {
      // Transport the nodes (returns false if 0 electrons in gap)
      if (!TransportTimeStep()) {
        break;
      }
    }
  }

  if (!m_vNElectronEvolution.empty()) {
    // determine maximal size of electron maxSize at time maxTime.
    auto maxSize = std::max_element(m_vNElectronEvolution.begin(),
                                    m_vNElectronEvolution.end(),
                                    [](const std::pair<double, long> &p1,
                                       const std::pair<double, long> &p2) {
                                      return p1.second < p2.second;
                                    });
    double maxTime =
        m_time0 +
        m_dt *
            ((double)std::distance(m_vNElectronEvolution.begin(), maxSize) + 1);

    std::cout << m_className
              << "::StartGridAvalanche: Avalanche maximum size of "
              << maxSize->second << " electrons reached at " << maxTime
              << " ns.\n";

    std::cout << m_className
              << "::StartGridAvalanche: Final avalanche size (produced "
                 "positive charge) = "
              << m_nTotPosIons << " ended at t = " << m_time << " ns.\n";
  }
}

void AvalancheGridSpaceCharge::ExportGrid(const std::string &filename) {
  std::ofstream exportElectrons(filename + "_electrons.csv");
  if (!exportElectrons.is_open()) throw Exception("Error opening e- file");
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      exportElectrons << m_grid[iz][ir].nE << " ";
    }
    exportElectrons << "\n";
  }
  exportElectrons.close();

  std::ofstream exportPosIon(filename + "_posion.csv");
  if (!exportPosIon.is_open()) throw Exception("Error opening p+ file");
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      exportPosIon << std::floor(m_grid[iz][ir].nP) << " ";
    }
    exportPosIon << "\n";
  }
  exportPosIon.close();

  std::ofstream exportNegIon(filename + "_negion.csv");
  if (!exportNegIon.is_open()) throw Exception("Error opening n- file");
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      exportNegIon << std::floor(m_grid[iz][ir].nN) << " ";
    }
    exportNegIon << "\n";
  }
  exportNegIon.close();

  std::ofstream exportZField(filename + "_eFieldZ.csv");
  if (!exportZField.is_open()) throw Exception("Error opening E_z file");
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      int gasGap = m_grid[iz][ir].gasGapIndex;
      const double ez = m_grid[iz][ir].ez + m_ezBkg[gasGap];
      exportZField << std::round(ez) << " ";
    }
    exportZField << "\n";
  }
  exportZField.close();

  std::ofstream exportRField(filename + "_eFieldR.csv");
  if (!exportRField.is_open()) throw Exception("Error opening E_r file");
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      const double er = m_grid[iz][ir].er;
      exportRField << std::round(er) << " ";
    }
    exportRField << "\n";
  }
  exportRField.close();

  std::ofstream exportMagField(filename + "_MagField.csv");
  if (!exportMagField.is_open()) throw Exception("Error opening E_r file");
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      int gasGap = m_grid[iz][ir].gasGapIndex;
      const double emag = Mag(m_grid[iz][ir].ez + m_ezBkg[gasGap], m_grid[iz][ir].er);
      exportMagField << std::round(emag) << " ";
    }
    exportMagField << "\n";
  }
  exportMagField.close();

  if (m_bDebug) {
    std::cout << m_className << "::ExportGrid: Grids exported.\n";
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Private Section:
bool AvalancheGridSpaceCharge::SnapTo2dGrid(const double x, const double y,
                                            const double z, const long n,
                                            const int gasLayer) {
  // Snap electron from AvalancheMicroscopic to the predefined grid
  if (m_grid.empty()) throw Exception("Grid is not defined");

  // y in micro is z in grid space-charge
  const double r = Mag(x - m_vCoNGasLayer[gasLayer][0], 
                       z - m_vCoNGasLayer[gasLayer][2]);
  int iZ = (int)std::round((y - m_zGrid.front()) * m_zInvStep);
  int iR = (int)std::round(r * m_rInvStep);

  if (m_bDebug) {
    std::cout << m_className << "::SnapTo2dGrid: iz = " << iZ << ", ir = " << iR
              << ".\n";
  }

  if (iZ < 0 || iZ > m_zSteps || iR < 0 || iR > m_rSteps) {
    if (m_bDebug) {
      std::cerr << m_className
                << "::SnapTo2dGrid: Point is outside the grid.\n";
    }
    return false;
  }

  // Add point to the grid.
  // When snapping the electron to the grid the distance traveled can yield
  // additional electrons or get attached. (depends on if against E field or
  // along ...). e-field is along y (micro)
  double step = m_zGrid[iZ] - y;
  // determine if against (ok) or with e field (not ok):
  int against = (step > 0 && m_ezBkg[gasLayer] < 0) ||
                (step < 0 && m_ezBkg[gasLayer] > 0);

  // sanity check
  if (m_grid[iZ][iR].gasGapIndex != gasLayer) {
    std::cerr << m_className
              << "::SnapTo2dGrid: Gas layer index does not match.\n";
    return false;
  }

  if (!against) {
    m_grid[iZ][iR].nE += n;
    m_nTotElectron += n;
    if (m_bDebug)
      std::cerr << m_className
                << "::SnapTo2dGrid: snap along e-field, continue.\n";
    return true;
  }

  // make step positive
  long nEOut;
  double nPOut, nNOut;
  GetAvalancheSizeFromStep(std::abs(step), n, m_grid[iZ][iR].townsend,
                           m_grid[iZ][iR].attachment, nEOut, nPOut, nNOut);
  if (nEOut == 0) {
    if (m_bDebug)
      std::cerr << m_className << "::SnapTo2dGrid: e- from " << n
                << " to 0 -> cancel.\n";
    return false;
  }

  m_grid[iZ][iR].nE += nEOut;
  m_grid[iZ][iR].nP += nPOut;
  m_grid[iZ][iR].nN += nNOut;
  m_nTotElectron += nEOut;
  m_nTotPosIons += (long)nPOut;

  if (m_bDebug) {
    std::cout << m_className << "::SnapTo2dGrid: e- from " << n << " to "
              << nEOut << " p+: " << nPOut << " n-: " << nNOut << ".\n"
              << "    Snapped to (z, r) = (" << y << " -> " << m_zGrid[iZ]
              << ", " << r << " -> " << m_rGrid[iR] << ").\n";
  }
  return true;
}

void AvalancheGridSpaceCharge::Prepare2dMesh() {
  // check if sensor is defined
  if (!m_sensor) {
    std::cerr << m_className
              << "::Prepare2dMesh: Sensor is not defined. Abort.\n";
  }

  // get a point (Y global coordinate) in each gas gap
  int n = m_vIndexGasGaps.size();
  m_ezBkg.resize(n);

  if (m_pp) {
    m_vYPointInGasGap.resize(n);
    for (int iz = 0; iz <= m_zSteps; iz++) {
      // Determine layer and gas gap.
      const int layerIndex = m_pp->GetLayer(m_zGrid[iz]);
      const int k = GetGasGapNumber(layerIndex);
      if (k != -1 && k < n) m_vYPointInGasGap[k] = m_zGrid[iz];
    }
  }

  std::vector<double> alpha(n), eta(n), vd(n), dSigmaL(n), dSigmaT(n), wv(n),
      wr(n), alphaPT(n), etaPT(n);
  double e[3];
  int status;
  Medium *m = nullptr;
  // iterate through the gas gaps
  for (int k = 0; k < n; k++) {
    m_sensor->ElectricField(0, m_vYPointInGasGap[k], 0, e[0], e[1], e[2], m,
                            status);

    if (status != 0) {
      std::cerr
          << m_className
          << "::Prepare2dMesh: Cannot estimate background field for gas gap "
          << k + 1 << ".\n";
    }

    // one expects (ComponentParallelPlate) that the electric field is pointing
    // along y-axis
    //  i.e. Z-axis in our coordinate system.
    m_ezBkg[k] = e[1];
    GetSwarmParameters(0., m_vYPointInGasGap[k], 0., std::abs(e[1]), 
                       alpha[k], eta[k], vd[k], dSigmaL[k], dSigmaT[k], 
                       wv[k], wr[k], alphaPT[k], etaPT[k]);

    // print-out to double-check the swarm parameters
    std::cout << m_className << "::Prepare2dMesh:\n"
              << "  Gas gap " << k + 1 << "\n"
              << "     Ez: " << m_ezBkg[k] << " (V/cm)\n"
              << "     alphaSST: " << alpha[k] << " (1/cm)\n"
              << "     alphaPT:  " << alphaPT[k] << " (1/cm)\n"
              << "     etaSST: " << eta[k] << " (1/cm)\n"
              << "     etaPT:  " << etaPT[k] << " (1/cm)\n"
              << "     drift velocity (Wv): " << vd[k] << " (cm/ns)\n"
              << "     Wr (!= Wv): " << wr[k] << " (cm/ns).\n";
  }

  // Set up mesh
  m_grid.resize(m_zSteps + 1);
  m_zGasGapBoundaries.resize(n);
  for (int iz = 0; iz <= m_zSteps; iz++) {
    m_grid[iz].resize(m_rSteps + 1);
    // Determine the gas gap.
    int k = 0;
    if (m_pp) {
      const int layerIndex = m_pp->GetLayer(m_zGrid[iz]);
      if (layerIndex >= 0) k = GetGasGapNumber(layerIndex);
    }
    if (k != -1) {
      // store z index for gap k
      m_zGasGapBoundaries[k].push_back(iz);
    }
    for (int ir = 0; ir <= m_rSteps; ir++) {
      // Set gas gap index.
      m_grid[iz][ir].gasGapIndex = k;
      // Continue if nodes are not in a gas gap.
      if (k == -1) continue;
      // Set swarm parameters.
      m_grid[iz][ir].townsend = alpha[k];
      m_grid[iz][ir].attachment = eta[k];
      m_grid[iz][ir].vd = vd[k];  //< magnitude! direction against E field
      m_grid[iz][ir].dSigmaL = dSigmaL[k];
      m_grid[iz][ir].dSigmaT = dSigmaT[k];
      m_grid[iz][ir].wv = wv[k];
      m_grid[iz][ir].wr = wr[k];
      m_grid[iz][ir].townsendPT = alphaPT[k];
      m_grid[iz][ir].attachmentPT = etaPT[k];
    }
  }

  // set anode in each gas gap
  for (int k = 0; k < n; k++) {
    int izMin = m_zGasGapBoundaries[k].front();
    int izMax = m_zGasGapBoundaries[k].back();
    int izAnode = (m_ezBkg[k] > 0) ? izMin : izMax;
    for (int ir = 0; ir <= m_rSteps; ir++) {
      m_grid[izAnode][ir].anode = true;
    }
  }

  // we define the time step as dz / max(Wr)
  m_dt = m_zStepSize / *std::max_element(wr.begin(), wr.end());

  if (m_bDebug) {
    std::cout << m_className << "::Prepare2dMesh: Time step per loop: " << m_dt
              << " ns.\n";
  }
}

void AvalancheGridSpaceCharge::PrepareElectronsFromMicroscopicAvalanche() {
  double tMicro = 0;
  long neTotal = 0;
  for (int k = 0; k < (int)m_vIndexGasGaps.size(); k++) {
    // calculate middle coord of electron cloud and add all electrons to the
    // grid/mesh per gas gap.
    auto neGasGap = m_vElectrons[k].size();

    // continue if no electron in the gas gap
    if (neGasGap <= 0) {
      // HS: not sure it's a good idea to initialize with NAN.
      m_vCoNGasLayer.push_back({NAN, NAN, NAN});
      continue;
    }
    neTotal += neGasGap;
    double xMicro = 0, yMicro = 0, zMicro = 0;

    // set time and center of electron number
    for (const auto &electron : m_vElectrons[k]) {
      xMicro += electron.x;
      yMicro += electron.y;
      zMicro += electron.z;
      tMicro += electron.t;
    }
    m_vCoNGasLayer.push_back({xMicro / (double)neGasGap,
                              yMicro / (double)neGasGap,
                              zMicro / (double)neGasGap});

    if (m_bDebug) {
      std::cout
          << m_className
          << "::PrepareElectronsFromMicroscopicAvalanche: mean center gas gap "
          << k + 1 << " X = (" << m_vCoNGasLayer[k][0] << ", "
          << m_vCoNGasLayer[k][1] << ", " << m_vCoNGasLayer[k][2] << ").\n";
    }
  }

  // Set the start time.
  m_time0 = tMicro / (double)neTotal;
  m_time = m_time0;

  // prepare the mesh (needs m_vCoN)
  Prepare2dMesh();

  // place all electrons onto the grid
  for (int k = 0; k < (int)m_vIndexGasGaps.size(); k++)
    for (auto &electron : m_vElectrons[k]) {  //< does not proceed if no
                                              // electrons in respective gap
      if (SnapTo2dGrid(electron.x, electron.y, electron.z, 1, k) && m_bDebug) {
        std::cout << m_className
                  << "::PrepareElectronsFromMicroscopicAvalanche: "
                  << "Electron added in gas gap " << k + 1 << "\n"
                  << "      at (x,y,z) =  (" << electron.x << "," << electron.y
                  << "," << electron.z << ").\n";
      }
    }
}

void AvalancheGridSpaceCharge::GetSwarmParameters(
    const double x, const double y, const double z, const double emag, 
    double &alpha, double &eta, double &vd, double &dSigmaL, double &dSigmaT, 
    double &wv, double &wr, double &alphaPT, double &etaPT) const {
  if (m_bDebug && false)
    std::cerr << m_className
              << "::GetSwarmParameters: Getting parameters at "
                 "|E| = "
              << emag << ".\n";

  // medium from sensor
  Medium *m = m_sensor->GetMedium(x, y, z);

  // alpha
  m->ElectronTownsend(0., emag, 0., 0., 0., 0., alpha);
  // eta
  m->ElectronAttachment(0., emag, 0., 0., 0., 0., eta);

  // mag of velocity
  double vx, vy, vz;
  m->ElectronVelocity(0., emag, 0., 0., 0., 0., vx, vy, vz);
  vd = std::sqrt(vx * vx + vy * vy + vz * vz);  //< Wv in Magboltz
  // wv, wr <- take only Wr and 'vd' for Wv
  wr = 0.;
  if (!m_bUseTOF ||
      !m->ElectronVelocityFluxBulk(0., emag, 0., 0., 0., 0., wv, wr) ||
      wr < Small) {
    wr = vd;
  }
  // rates, if not available we take (alpha-eta)SST and ratio of alpha/eta =
  // Rion/Ratt Rion-Ratt = Reff (tagashira eq.)
  //  -> Rion converged to a rate with Wr and DL (using the one equation),
  //  alphaSST is either from SST and if not converged from magboltz itself.
  double rion = 0, ratt = 0;
  if (!m_bUseTOF ||
      !m->ElectronTOFIonisation(0., emag, 0., 0., 0., 0., rion) ||
      !m->ElectronTOFAttachment(0., emag, 0., 0., 0., 0., ratt)) {
    if (m_bDebug) {
      std::cerr << m_className
                << "::GetSwarmParameters: TOF Rates not available.\n";
    }

    // Diffusionless approximation
    rion = alpha * wr;
    ratt = eta * wr;
  }
  // calculate alpha/eta PT
  alphaPT = rion / wr;
  etaPT = ratt / wr;

  // diffusion coefficients
  m->ElectronDiffusion(0., emag, 0., 0., 0., 0., dSigmaL, dSigmaT);

  // print (and information about units!)
  if (m_bDebug && false) {
    std::cout << m_className << "::GetSwarmParameters:\n"
              << "  Townsend = " << alpha << " [1/cm], Attachment = " << eta
              << " [1/cm], Flux Velocity = "
              << vd
              //                      << " [cm/ns], Wv = " << wv
              << " [cm/ns], Bulk velocity = " << wr << " [cm/ns].\n"
              << "  TOF Ionization rate = " << alphaPT
              << " [1/ns], TOF Attachment rate = " << etaPT << "[1/ns].\n"
              << "  Longitudinal Diffusion = " << dSigmaL << " [sqrt(cm)],"
              << " Transversal Diffusion = " << dSigmaT << " [sqrt(cm)].\n";
  }
}

bool AvalancheGridSpaceCharge::TransportTimeStep() {
  // Propagate grid nodes by one time step with updated electric fields 
  // (Lippmann et al. approach).
  if (m_nTotElectron <= 0) return false;

  if (m_bDebug) {
    std::cout << m_className << "::TransportTimeStep: Start time: " << m_time
              << "\n";
  }

  if (m_bSpaceCharge && m_nTotElectron > 1e5) {
    // Clear existing rings.
    for (auto & ringsystem : m_vRingSystems) {
      ringsystem.ClearActiveRings();
      ringsystem.UpdateCentre(0., 0.);
    }
    
    for (int iz = 0; iz <= m_zSteps; iz++) {
      for (int ir = 0; ir <= m_rSteps; ir++) {
        const double q = -m_grid[iz][ir].nE + m_grid[iz][ir].nP - m_grid[iz][ir].nN;
        // If there is enough charge, count as an active node
        if (std::abs(q) < 0.1) continue;
        double zf = m_zGrid[iz];
        double rf = m_rGrid[ir];

        int gasGapIndex = m_grid[iz][ir].gasGapIndex;
        // Add the ring to the correct system: need the index of the gasgap
        // Direct charge interaction
        m_vRingSystems[gasGapIndex].AddChargedRing(rf, zf, 0., q); 

        if (m_fieldOption == FieldOption::Mirror) {
          // Assume symmetric single-layer RPC with resistive layers of
          // equal permittivity.
          if (m_vIndexGasGaps.size() > 1) {
            throw std::runtime_error(
                "::TransportTimeStep: Mirror charge option implemented but not tested for "
                "MRPC.");
          }

          // HS: this can be done at initialization time...
          // get epsilon value from neighboring layer (assume both layers have same
          // eps)
          int IndexOfRightLayer = m_vIndexGasGaps[gasGapIndex] + 1;
          // int IndexOfLeftLayer = m_vIndexGasGaps[gasGapIndex] - 1;
          double eps = m_pp->GetPermittivityFromLayer(IndexOfRightLayer);
          double alpha12 = (1. - eps) / (1. + eps);

          // Obtain bounds of current gas gap
          double zTop, zBottom;
          m_pp->getZBoundFromLayer(m_vIndexGasGaps[gasGapIndex], zTop, zBottom);

          // mirror charge interaction
          for (int i = 0; i < m_iFieldApprox; i++) {
            if (i == 0) {
              // 2a, alpha12 = delta_Q
              double zf0 = zf + 2. * (zTop - zf);
              m_vRingSystems[gasGapIndex].AddChargedRing(rf, zf0, 0., q * alpha12);
  
              // -2a', alpha12 = delta_Q
              zf0 = zf + 2. * (zBottom - zf);
              m_vRingSystems[gasGapIndex].AddChargedRing(rf, zf0, 0., q * alpha12);
            } else if (i == 1) {
              // TODO: higher order mirror charges
            } else {
              continue;
            }
          }
        }
      }
    }
  } 

  // choose MC or Mean version depending on m_bMC; total electron > 1e5
  std::function<void(double, const long, const double, const double, long &,
                     double &, double &)>
      AvalancheGain = GetAvalancheSizeFromStep;
  if (!m_bMC && m_nTotElectron > 1e5) {
    AvalancheGain = GetMeanAvalancheSizeFromStep;
    // TODO: Diffusion
  }

  // update the nodes for the next run (SC-field and swarm parameters)
  // MRPC: SC-effect only within each gas gap and option="coulomb"
  if (m_bSpaceCharge && m_nTotElectron > 1e5) {
    double dummy;  // dummy E field component as we are using a 2d grid
    Medium *m = nullptr;
    int stat;

    for (int iz = 0; iz <= m_zSteps; iz++) {
      double zi = m_zGrid[iz];
      // continue if not in gas gap
      int gasGap = m_grid[iz][0].gasGapIndex;
      if (gasGap == -1) continue;
      for (int ir = 0; ir <= m_rSteps; ir++) {
        double ri = m_rGrid[ir];
        auto &nd = m_grid[iz][ir];

        // reset local fields at node
        nd.ez = 0;
        nd.er = 0;

        // Skip if there are no electrons or if we are at an anode.
        if (nd.nE < 1 || nd.anode) continue;

        // update space charge field on each node
        int gasGapIndex = m_grid[iz][ir].gasGapIndex;
        m_vRingSystems[gasGapIndex].ElectricField(ri, zi, 0., nd.er, nd.ez, dummy, m, stat);
        
        // check if local field reaches background field values.
        const double emag = Mag(nd.ez + m_ezBkg[gasGap], nd.er);
        if (emag - std::abs(m_ezBkg[gasGap]) >=
                m_fStreamerK * std::abs(m_ezBkg[gasGap]) &&
            !m_bFieldK) {
          std::cout << m_className << ":TransportTimeStep:\n"
                    << "    Space-charge field reached "
                    << std::to_string(int(m_fStreamerK * 100))
                    << "% of background field in gas gap " << gasGap + 1
                    << "\n";
          // TODO: Total electrons in gas gap "gasGap"
          m_lElectronsK = m_vNElectronEvolution.back().second;
          m_bFieldK = true;
        }

        // calculate the swarm parameters
        GetSwarmParameters(0., zi, 0., emag, nd.townsend, nd.attachment, nd.vd,
                           nd.dSigmaL, nd.dSigmaT, nd.wv, nd.wr, nd.townsendPT,
                           nd.attachmentPT);

        // get new step distance
        double step = std::abs(nd.wr * m_dt);

        // adaptive time stepping (this routine takes the smallest dt needed for
        // the current sc-field)
        if (m_bAdaptiveTime && step >= 2. * m_zStepSize) {
          // reset dt
          double dtPrev = m_dt;
          m_dt = m_zStepSize / nd.wr;

          if (m_bDebug) {
            std::cout << m_className << "::TransportTimeStep: Changed dt from "
                      << dtPrev << " to: " << m_dt << "\n"
                      << "      due to step size: " << step
                      << " bulk velocity: " << nd.wr << "\n"
                      << "      electric field: " << emag
                      << " alpha: " << nd.townsendPT
                      << " eta: " << nd.attachmentPT << "\n"
                      << "      diffusion longitudinal/transversal: "
                      << nd.dSigmaL << " " << nd.dSigmaT << "\n";
            ExportGrid("TIME_STEP_ADAPTION_" + std::to_string(m_dt));
          }
        }
      }
    }
  }

  // finish if stop reached and set at 100 * K %
  if (m_bStopAtK && m_bFieldK) return false;

  // transport the electrons with the update sc-field, swarm parameter and dt
  for (int iz = 0; iz <= m_zSteps; iz++) {
    // continue if not in gas gap
    const int gasGap = m_grid[iz][0].gasGapIndex;
    if (gasGap == -1) continue;
    bool saturated = false;
    if (!m_bSpaceCharge &&
        (std::find(m_vSaturatedGaps.begin(), m_vSaturatedGaps.end(),
                   gasGap) != m_vSaturatedGaps.end())) {
      saturated = true;
    }
    for (int ir = 0; ir <= m_rSteps; ir++) {
      auto &nd = m_grid[iz][ir];

      // Skip if we are at the anode or if there are no electrons.
      if (nd.anode || nd.nE < 1) continue;

      // update step distance
      double step = std::abs(nd.wr * m_dt);

      // calculate new avalanche size at X + step
      // HS: use a vector<bool> to keep track of which gaps are saturated?
      long nEOut = 0;
      double nPOut = 0.;
      double nNOut = 0.;
      if (saturated) {
        // Saturated case, don't evolve electrons in size
        nEOut = nd.nE;
        nPOut = 0,
        nNOut = 0;  //< strictly this is completely wrong because
                    // SC-bremsung creates huge amounts of ions
      } else {
        AvalancheGain(step, nd.nE, nd.townsendPT, nd.attachmentPT,
                      nEOut, nPOut, nNOut);
      }
      m_nTotPosIons += std::round(nPOut);

      // calculate steps against electric field i.e. correct sign.
      const double emag = Mag(nd.ez + m_ezBkg[gasGap], nd.er);
      double stepZ = step * (-(nd.ez + m_ezBkg[gasGap]) / emag);
      const double stepR = step * (-(nd.er) / emag);

      if (m_bDiffusion) {
        // correct the stepping from diffusion + charge distribution
        DiffuseTimeStep(step, emag, nEOut, std::round(nPOut), std::round(nNOut),
                        iz, ir, gasGap);
      } else {
        // calculate steps and distribute charges (no diffusion)
        DistributeCharges(nEOut, std::round(nPOut), std::round(nNOut), 
                          iz, ir, stepZ, stepR, gasGap);
      }

      if (m_sensor->GetNumberOfElectrodes() > 0) {
        // For now, use a single charge at phi = 0.
        // TODO: 
        // - Discretize phi and add signal for each slice with nElectron / M.
        // - At anode, the signal from diffusion is due to bounded plane not
        //   net 0 because diffusion tends more backwards, since forward they 
        //   reach the boundary at earlier distance.
        constexpr double cphi = 1.;
        constexpr double sphi = 0.;
        const double r0 = m_rGrid[ir];
        const double x0 = m_vCoNGasLayer[gasGap][0] + r0 * cphi;
        const double y0 = m_zGrid[iz];
        const double z0 = m_vCoNGasLayer[gasGap][2] - r0 * sphi;

        // z-step outside gasGap domain, resize to stepZ = Anode - Current
        int izMin = m_zGasGapBoundaries[gasGap].front();
        int izMax = m_zGasGapBoundaries[gasGap].back();
        if (m_zGrid[iz] + stepZ < m_zGrid[izMin]) {
          stepZ = (m_zGrid[izMin] - m_zGrid[iz]);
        } else if (m_zGrid[iz] + stepZ > m_zGrid[izMax]) {
          stepZ = (m_zGrid[izMax] - m_zGrid[iz]);
        }

        const double r1 = m_rGrid[ir] + stepR;
        const double x1 = m_vCoNGasLayer[gasGap][0] + r1 * cphi;
        const double y1 = m_zGrid[iz] + stepZ;
        const double z1 = m_vCoNGasLayer[gasGap][2] - r1 * sphi;

        // Induced current from flux drift velocity i.e. introduce weight factor
        //< 1 if (Wv = velocity): Wr = flux
        const double weight = nd.vd / nd.wr;  
        m_sensor->AddSignalWeightingPotential(
            -weight, {m_time, m_time + m_dt}, {{x0, y0, z0}, {x1, y1, z1}},
            {(double)nd.nE, (double)nEOut});
      }
    }
  }

  // propagate Grid time (after above for loops due to the adaptive time
  // stepping)
  m_time += m_dt;

  // sums electrons left in gap (count per gap)
  std::vector<long> eOnGrid(m_vIndexGasGaps.size(), 0);

  // update nodes with transported electrons
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      // Get node
      auto &nd = m_grid[iz][ir];
      int gasGap = nd.gasGapIndex;
      if (nd.anode && m_bStick) {
        // sticky anode: electron stay and holder electrons add it up
        nd.nE += nd.nEHolder;
      } else {
        // update node with electrons from holder
        nd.nE = nd.nEHolder;
      }
      // ions add up, also at the anode
      nd.nP += nd.nPHolder;
      nd.nN += nd.nNHolder;

      // reset node Holder
      nd.nEHolder = 0;
      nd.nPHolder = 0;
      nd.nNHolder = 0;

      // add electrons if they are not stuck
      if (!(nd.anode && m_bStick)) eOnGrid[gasGap] += nd.nE;
    }
  }
  // add total electrons in gap to grid and to evolution vector
  m_nTotElectron = std::accumulate(eOnGrid.begin(), eOnGrid.end(), 0.);
  m_vNElectronEvolution.push_back(std::make_pair(m_time, m_nTotElectron));

  // determine saturated gaps at each time step
  // clear: anode-absorption activates avalanche to grow again
  m_vSaturatedGaps.resize(0);
  for (int k = 0; k < (int)eOnGrid.size(); k++) {
    if (!m_bSpaceCharge && eOnGrid[k] > m_lNCrit) {
      m_vSaturatedGaps.push_back(k);
    }
    if (m_bDebug) {
      std::cout << m_className
                << "::TransportTimeStep: Electrons active on grid in gas gap "
                << k + 1 << ": " << eOnGrid[k] << "\n";
    }
  }

  return true;
}

void AvalancheGridSpaceCharge::DiffuseTimeStep(
    const double dx, const double emag,
    const long nE, const double nP, const double nN,
    const int iz, const int ir, const int gap) {
  // Add diffusion onto the step dx

  // Minimum number of groups (same values as Lippmann).
  constexpr long nMinGroups = 50; 
  // Same values as Lippmann & Riegler
  constexpr std::array<long, 10> groupSizes = {
    1500, 800, 400, 200, 100,
      50,   20,  10,  5,   2};  

  // Split in at least nMinGroups groups of size groupSize.
  long rest = 0, groups = 0, groupSize = 0;
  for (const auto size : groupSizes) {
    if (nE > nMinGroups * size) {
      groupSize = size;                     //< size of a single group
      groups = std::floor(nE / groupSize);  //< real # groups
      rest = nE - groups * groupSize;       //< rest
      break;
    }
  }

  // If there are too few electrons we diffuse each electron individually.
  if (nE <= nMinGroups * groupSizes.back()) {
    groupSize = 1;
    groups = nE;
    rest = 0;
  }

  // Calculate diffusion and add to transport step.
  double sinTheta = 0.;
  double cosTheta = 1.;
  auto &nd = m_grid[iz][ir];
  if (emag > 1.e-8) {
    const double einv = 1. / emag;
    cosTheta = (-(nd.ez + m_ezBkg[gap]) * einv);
    sinTheta = (-(nd.er) * einv);
  }

  double f = (double)groupSize / (double)nE;
  const double r = m_rGrid[ir];
  const double sqrtdx = std::sqrt(dx);

  for (int group = 0; group < groups; group++) {
    // In the last loop we add the rest to the groupSize.
    if (group == groups - 1) {
      groupSize += rest;
      f = (double)groupSize / (double)nE;
    }
    // Diffuse each group as if it is a particle.
    // Calculate a diffusion step in a local coordinate system (U,V,W)
    // where W is along the E field, V is along e_phi,
    // and U is perpendicular to V and W.
    const double dU = RndmGaussian(0, nd.dSigmaT * sqrtdx);
    const double dV = RndmGaussian(0, nd.dSigmaT * sqrtdx);
    const double dW = RndmGaussian(dx, nd.dSigmaL * sqrtdx);

    // Transform to avalanche coordinate system
    // (Z,R,Y) where R mimics an X axis and Y is perpendicular to R and Z
    const double dX = cosTheta * dU + sinTheta * dW;
    // dY = dV
    // Sign seems correct due to sign in cos- and sinTheta
    const double dZ = cosTheta * dW - sinTheta * dU;  

    // calculate the change of radius
    // sign correct and stepR >= -r
    const double stepR = std::sqrt((r + dX) * (r + dX) + dV * dV) - r; 

    // Distribute nodes and add electrons/ions to Holder.
    // Fractional ion number is allowed otherwise loss of ions.
    DistributeCharges(groupSize, nP * f, nN * f, iz, ir, dZ, stepR, gap);
  }
}

void AvalancheGridSpaceCharge::DistributeCharges(long nE, double nP,
                                                 double nN, int iz, int ir,
                                                 double stepZ, double stepR,
                                                 int gasGap) {
  // distributes the charges from a movement in Z and R direction
  // Compute the effective step size in r-direction
  stepR = std::abs(m_rGrid[ir] + stepR) - m_rGrid[ir];
  double zRatio = (m_zGrid[iz] + stepZ - m_zGrid[0]) * m_zInvStep;
  double rRatio = std::abs(m_rGrid[ir] + stepR) * m_rInvStep;  // can travel through r=0

  int iz1, signZstep;
  if (stepZ < 0) {
    iz1 = (int)ceil(zRatio);
    signZstep = -1;
  } else if (stepZ == 0) {
    iz1 = iz;
    signZstep = 0;
  } else {
    iz1 = (int)floor(zRatio);
    signZstep = +1;
  }

  int ir1, signRstep;
  if (stepR < 0) {
    ir1 = (int)ceil(rRatio);
    if (ir1 == 0) {
      // otherwise it travels to ir = -1
      signRstep = +1;
    } else {
      signRstep = -1;
    }
  } else if (stepR == 0) {
    // no movement in r direction
    ir1 = ir;
    signRstep = 0;
  } else {
    ir1 = (int)floor(rRatio);
    signRstep = +1;
  }

  // 4 point approximation:
  int iz2 = iz1 + signZstep;
  // always: ir2 >= 0
  int ir2 = ir1 + signRstep;

  const double bz = std::abs(zRatio - (double)iz1);
  const double az = 1. - bz;
  const double br = std::abs(rRatio - (double)ir1);
  const double ar = 1. - br;

  if (az < 0 || az > 1) throw std::runtime_error("az not in range");
  if (bz < 0 || bz > 1) throw std::runtime_error("bz not in range");
  if (ar < 0 || ar > 1) throw std::runtime_error("ar not in range");
  if (br < 0 || br > 1) throw std::runtime_error("br not in range");

  // check if still in grid else place at boundaries (will be absorbed in next
  // step)?
  int izMin = m_zGasGapBoundaries[gasGap].front();
  int izMax = m_zGasGapBoundaries[gasGap].back();
  if (iz1 < izMin) iz1 = izMin;
  if (iz2 < izMin) iz2 = izMin;
  if (iz1 > izMax) iz1 = izMax;
  if (iz2 > izMax) iz2 = izMax;

  if (ir1 > m_rSteps) ir1 = m_rSteps;
  if (ir2 > m_rSteps) ir2 = m_rSteps;

  // add to the nodes the electrons travelled to (into nEHolder as they
  // will mix)
  if (nE > 200) {
    m_grid[iz1][ir1].nEHolder += (long)std::round(nE * az * ar);
    m_grid[iz1][ir2].nEHolder += (long)std::round(nE * az * br);
    m_grid[iz2][ir1].nEHolder += (long)std::round(nE * bz * ar);
    m_grid[iz2][ir2].nEHolder += (long)std::round(nE * bz * br);

    // add positive ions to the nodes (smeared values allowed)
    m_grid[iz1][ir1].nPHolder += nP * az * ar;
    m_grid[iz1][ir2].nPHolder += nP * az * br;
    m_grid[iz2][ir1].nPHolder += nP * bz * ar;
    m_grid[iz2][ir2].nPHolder += nP * bz * br;

    // add negative ions to the nodes (smeared values allowed)
    m_grid[iz1][ir1].nNHolder += nN * az * ar;
    m_grid[iz1][ir2].nNHolder += nN * az * br;
    m_grid[iz2][ir1].nNHolder += nN * bz * ar;
    m_grid[iz2][ir2].nNHolder += nN * bz * br;

  } else {
    // Too few electrons.
    // Only move them to 1 node (instead of 4).
    iz1 = (az >= bz) ? iz1 : iz2;
    ir1 = (ar >= br) ? ir1 : ir2;

    m_grid[iz1][ir1].nEHolder += nE;
    m_grid[iz1][ir1].nPHolder += nP;
    m_grid[iz1][ir1].nNHolder += nN;
  }
}

double AvalancheGridSpaceCharge::GetMeanDistance() {
  // Returns mean distance of electrons on the whole grid (doesn't work for
  // MRPCs)
  long nofElectrons = 0;
  double z = 0.;
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      const auto ne = m_grid[iz][ir].nE;
      if (ne < 0.5) continue;
      nofElectrons += ne;
      z += m_zGrid[iz] * ne;
    }
  }
  return z / (double)nofElectrons;
}

void AvalancheGridSpaceCharge::SetRingSystems() {
  // We add a charged ring system for each gas gap
  if (!m_pp)
    std::cerr << m_className
              << "::SetRingSystems: Parallel plate improperly defined.\n";
  m_pp->IndexOfGasGaps(m_vIndexGasGaps);
  size_t n_gas_gaps = m_vIndexGasGaps.size();

  double horizontal_max = m_rGrid.back();
  double horizontal_min = -1. * horizontal_max;
  double vertical_min = m_zGrid.front();
  double vertical_max = m_zGrid.back();

  for (size_t i = 0; i < n_gas_gaps; ++i) {
    m_vRingSystems.emplace_back();
    const int layer_index = m_vIndexGasGaps[i];
    double y_bottom, y_top;
    m_pp->getZBoundFromLayer(layer_index, y_bottom, y_top);
    Medium *m = m_pp->GetMedium(0, 0.5 * (y_top - y_bottom) + y_bottom, 0);
    m_vRingSystems[i].SetMedium(m);
    m_vRingSystems[i].SetArea(horizontal_min, vertical_min, horizontal_min,
                              horizontal_max, vertical_max, horizontal_max);
    if (m_bDebug) m_vRingSystems[i].EnableDebugging();
  }
}

}  // namespace Garfield
