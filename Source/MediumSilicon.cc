#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Random.hh"
#include "Garfield/Utilities.hh"

namespace Garfield {

MediumSilicon::MediumSilicon()
    : Medium() {
  m_className = "MediumSilicon";
  m_name = "Si";

  SetTemperature(300.);
  SetDielectricConstant(11.9);
  Medium::SetAtomicNumber(14.);
  Medium::SetAtomicWeight(28.0855);
  Medium::SetMassDensity(2.329);

  m_driftable = true;
  m_ionisable = true;
  m_microscopic = true;

  m_w = 3.6;
  m_fano = 0.11;
 
  m_cb[0].eMin = 0.;
  m_cb[1].eMin = 1.05;
  m_cb[2].eMin = 2.24;
  m_cb[0].eFinal =  4.;
  m_cb[1].eFinal =  4.;
  m_cb[2].eFinal = 10.;
  for (size_t i = 0; i < 3; ++i) {
    m_cb[i].eStep = m_cb[i].eFinal / m_cb[i].nEnergySteps;
    m_cb[i].iMin = int(m_cb[i].eMin / m_cb[i].eStep) + 1;
  }
  // Effective masses.
  m_cb[0].mL = 0.916;
  m_cb[0].mT = 0.191;
  m_cb[1].mL = 1.59;
  m_cb[1].mT = 0.12;
  // Conduction effective masses.
  m_cb[0].mC = 3. / (1. / m_cb[0].mL + 2. / m_cb[0].mT);
  m_cb[1].mC = 3. / (1. / m_cb[1].mL + 2. / m_cb[1].mT);
  // Non-parabolicity parameters [1/eV].
  m_cb[0].alpha = 0.5;
  m_cb[1].alpha = 0.5;
  m_cb[0].nValleys = 6;
  m_cb[1].nValleys = 8;
  m_cb[2].nValleys = 1;
  for (size_t i = 0; i < 3; ++i) {
    const size_t nV = m_cb[i].nValleys;
    for (size_t j = 0; j < nV; ++j) m_cbIndex.push_back(i);
  }

  m_vb.eFinal =  8.5;
  m_vb.eStep = m_vb.eFinal / m_vb.nEnergySteps;

  // Load the density of states table.
  InitialiseDensityOfStates();
}

void MediumSilicon::SetDoping(const char type, const double c) {
  if (toupper(type) == 'N') {
    m_dopingType = 'n';
    if (c > Small) {
      m_dopingConcentration = c;
    } else {
      std::cerr << m_className << "::SetDoping:\n"
                << "    Doping concentration must be greater than zero.\n"
                << "    Using default value for n-type silicon "
                << "(10^12 cm-3) instead.\n";
      m_dopingConcentration = 1.e12;
    }
  } else if (toupper(type) == 'P') {
    m_dopingType = 'p';
    if (c > Small) {
      m_dopingConcentration = c;
    } else {
      std::cerr << m_className << "::SetDoping:\n"
                << "    Doping concentration must be greater than zero.\n"
                << "    Using default value for p-type silicon "
                << "(10^18 cm-3) instead.\n";
      m_dopingConcentration = 1.e18;
    }
  } else if (toupper(type) == 'I') {
    m_dopingType = 'i';
    m_dopingConcentration = 0.;
  } else {
    std::cerr << m_className << "::SetDoping:\n"
              << "    Unknown dopant type (" << type << ").\n"
              << "    Available types are n, p and i (intrinsic).\n";
    return;
  }

  m_isChanged = true;
}

void MediumSilicon::GetDoping(char& type, double& c) const {
  type = m_dopingType;
  c = m_dopingConcentration;
}

void MediumSilicon::SetTrapCrossSection(const double ecs, const double hcs) {
  if (ecs < 0.) {
    std::cerr << m_className << "::SetTrapCrossSection:\n"
              << "    Capture cross-section [cm2] must non-negative.\n";
  } else {
    m_eTrapCs = ecs;
  }

  if (hcs < 0.) {
    std::cerr << m_className << "::SetTrapCrossSection:\n"
              << "    Capture cross-section [cm2] must be non-negative.n";
  } else {
    m_hTrapCs = hcs;
  }

  m_trappingModel = 0;
  m_isChanged = true;
}

void MediumSilicon::SetTrapDensity(const double n) {
  if (n < 0.) {
    std::cerr << m_className << "::SetTrapDensity:\n"
              << "    Trap density [cm-3] must be non-negative.\n";
  } else {
    m_eTrapDensity = n;
    m_hTrapDensity = n;
  }

  m_trappingModel = 0;
  m_isChanged = true;
}

void MediumSilicon::SetTrappingTime(const double etau, const double htau) {
  if (etau <= 0.) {
    std::cerr << m_className << "::SetTrappingTime:\n"
              << "    Trapping time [ns-1] must be positive.\n";
  } else {
    m_eTrapTime = etau;
  }

  if (htau <= 0.) {
    std::cerr << m_className << "::SetTrappingTime:\n"
              << "    Trapping time [ns-1] must be positive.\n";
  } else {
    m_hTrapTime = htau;
  }

  m_trappingModel = 1;
  m_isChanged = true;
}

bool MediumSilicon::ElectronVelocity(const double ex, const double ey,
                                     const double ez, const double bx,
                                     const double by, const double bz,
                                     double& vx, double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::ElectronVelocity:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_eVelE.empty()) {
    // Interpolation in user table.
    return Medium::ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }

  // Calculate the mobility.
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  const double mu = -ElectronMobility(emag);

  if (fabs(bx) < Small && fabs(by) < Small && fabs(bz) < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    Langevin(ex, ey, ez, bx, by, bz, mu, m_eHallFactor * mu, vx, vy, vz);
  }
  return true;
}

bool MediumSilicon::ElectronTownsend(const double ex, const double ey,
                                     const double ez, const double bx,
                                     const double by, const double bz,
                                     double& alpha) {
  alpha = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::ElectronTownsend:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_eAlp.empty()) {
    // Interpolation in user table.
    return Medium::ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  }

  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  alpha = ElectronAlpha(emag);
  return true;
}

bool MediumSilicon::ElectronAttachment(const double ex, const double ey,
                                       const double ez, const double bx,
                                       const double by, const double bz,
                                       double& eta) {
  eta = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::ElectronAttachment:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_eAtt.empty()) {
    // Interpolation in user table.
    return Medium::ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
  }

  switch (m_trappingModel) {
    case 0:
      eta = m_eTrapCs * m_eTrapDensity;
      break;
    case 1:
      double vx, vy, vz;
      ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
      eta = m_eTrapTime * sqrt(vx * vx + vy * vy + vz * vz);
      if (eta > 0.) eta = -1. / eta;
      break;
    default:
      std::cerr << m_className << "::ElectronAttachment: Unknown model. Bug!\n";
      return false;
  }

  return true;
}

bool MediumSilicon::HoleVelocity(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz, double& vx,
                                 double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::HoleVelocity:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_hVelE.empty()) {
    // Interpolation in user table.
    return Medium::HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }

  // Calculate the mobility.
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  const double mu = HoleMobility(emag);

  if (fabs(bx) < Small && fabs(by) < Small && fabs(bz) < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    Langevin(ex, ey, ez, bx, by, bz, mu, m_hHallFactor * mu, vx, vy, vz);
  }
  return true;
}

bool MediumSilicon::HoleTownsend(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz,
                                 double& alpha) {
  alpha = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::HoleTownsend:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_hAlp.empty()) {
    // Interpolation in user table.
    return Medium::HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }

  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  alpha = HoleAlpha(emag);
  return true;
}

bool MediumSilicon::HoleAttachment(const double ex, const double ey,
                                   const double ez, const double bx,
                                   const double by, const double bz,
                                   double& eta) {
  eta = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::HoleAttachment:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_hAtt.empty()) {
    // Interpolation in user table.
    return Medium::HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }

  switch (m_trappingModel) {
    case 0:
      eta = m_hTrapCs * m_hTrapDensity;
      break;
    case 1:
      double vx, vy, vz;
      HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
      eta = m_hTrapTime * sqrt(vx * vx + vy * vy + vz * vz);
      if (eta > 0.) eta = -1. / eta;
      break;
    default:
      std::cerr << m_className << "::HoleAttachment: Unknown model. Bug!\n";
      return false;
  }

  return true;
}

void MediumSilicon::SetLowFieldMobility(const double mue, const double muh) {
  if (mue <= 0. || muh <= 0.) {
    std::cerr << m_className << "::SetLowFieldMobility:\n"
              << "    Mobility must be greater than zero.\n";
    return;
  }

  m_eMobility = mue;
  m_hMobility = muh;
  m_hasUserMobility = true;
  m_isChanged = true;
}

void MediumSilicon::SetLatticeMobilityModelMinimos() {
  m_latticeMobilityModel = LatticeMobility::Minimos;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetLatticeMobilityModelSentaurus() {
  m_latticeMobilityModel = LatticeMobility::Sentaurus;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetLatticeMobilityModelReggiani() {
  m_latticeMobilityModel = LatticeMobility::Reggiani;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetDopingMobilityModelMinimos() {
  m_dopingMobilityModel = DopingMobility::Minimos;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetDopingMobilityModelMasetti() {
  m_dopingMobilityModel = DopingMobility::Masetti;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetSaturationVelocity(const double vsate,
                                          const double vsath) {
  if (vsate <= 0. || vsath <= 0.) {
    std::cout << m_className << "::SetSaturationVelocity:\n"
              << "    Restoring default values.\n";
    m_hasUserSaturationVelocity = false;
  } else {
    m_eSatVel = vsate;
    m_hSatVel = vsath;
    m_hasUserSaturationVelocity = true;
  }

  m_isChanged = true;
}

void MediumSilicon::SetSaturationVelocityModelMinimos() {
  m_saturationVelocityModel = SaturationVelocity::Minimos;
  m_hasUserSaturationVelocity = false;
  m_isChanged = true;
}

void MediumSilicon::SetSaturationVelocityModelCanali() {
  m_saturationVelocityModel = SaturationVelocity::Canali;
  m_hasUserSaturationVelocity = false;
  m_isChanged = true;
}

void MediumSilicon::SetSaturationVelocityModelReggiani() {
  m_saturationVelocityModel = SaturationVelocity::Reggiani;
  m_hasUserSaturationVelocity = false;
  m_isChanged = true;
}

void MediumSilicon::SetHighFieldMobilityModelMinimos() {
  m_highFieldMobilityModel = HighFieldMobility::Minimos;
  m_isChanged = true;
}

void MediumSilicon::SetHighFieldMobilityModelCanali() {
  m_highFieldMobilityModel = HighFieldMobility::Canali;
  m_isChanged = true;
}

void MediumSilicon::SetHighFieldMobilityModelReggiani() {
  m_highFieldMobilityModel = HighFieldMobility::Reggiani;
  m_isChanged = true;
}

void MediumSilicon::SetHighFieldMobilityModelConstant() {
  m_highFieldMobilityModel = HighFieldMobility::Constant;
}

void MediumSilicon::SetImpactIonisationModelVanOverstraetenDeMan() {
  m_impactIonisationModel = ImpactIonisation::VanOverstraeten;
  m_isChanged = true;
}

void MediumSilicon::SetImpactIonisationModelGrant() {
  m_impactIonisationModel = ImpactIonisation::Grant;
  m_isChanged = true;
}

void MediumSilicon::SetImpactIonisationModelMassey() {
  m_impactIonisationModel = ImpactIonisation::Massey;
  m_isChanged = true;
}

void MediumSilicon::SetImpactIonisationModelOkutoCrowell() {
  m_impactIonisationModel = ImpactIonisation::Okuto;
  m_isChanged = true;
}

bool MediumSilicon::SetMaxElectronEnergy(const double e) {
  if (e <= m_cb[2].eMin + Small) {
    std::cerr << m_className << "::SetMaxElectronEnergy:\n"
              << "    Requested upper electron energy limit (" << e
              << " eV) is too small.\n";
    return false;
  }

  m_cb[2].eFinal = e;
  // Determine the energy interval size.
  m_cb[2].eStep = m_cb[2].eFinal / m_cb[2].nEnergySteps;

  m_isChanged = true;

  return true;
}

double MediumSilicon::GetElectronEnergy(const double px, const double py,
                                        const double pz, double& vx, double& vy,
                                        double& vz, const int band) {
  // Effective masses
  double mx = ElectronMass, my = ElectronMass, mz = ElectronMass;
  // Energy offset
  double e0 = 0.;
  if (band >= 0 && band < m_cb[0].nValleys) {
    // X valley
    if (m_anisotropic) {
      switch (band) {
        case 0:
        case 1:
          // X 100, -100
          mx *= m_cb[0].mL;
          my *= m_cb[0].mT;
          mz *= m_cb[0].mT;
          break;
        case 2:
        case 3:
          // X 010, 0-10
          mx *= m_cb[0].mT;
          my *= m_cb[0].mL;
          mz *= m_cb[0].mT;
          break;
        case 4:
        case 5:
          // X 001, 00-1
          mx *= m_cb[0].mT;
          my *= m_cb[0].mT;
          mz *= m_cb[0].mL;
          break;
        default:
          std::cerr << m_className << "::GetElectronEnergy:\n"
                    << "    Unexpected band index " << band << "!\n";
          break;
      }
    } else {
      // Use the conduction effective mass.
      mx *= m_cb[0].mC;
      my *= m_cb[0].mC;
      mz *= m_cb[0].mC;
    }
  } else if (band < m_cb[0].nValleys + m_cb[1].nValleys) {
    // L valley, isotropic approximation
    e0 = m_cb[1].eMin;
    mx *= m_cb[1].mC;
    my *= m_cb[1].mC;
    mz *= m_cb[1].mC;
  } else if (band == m_cb[0].nValleys + m_cb[1].nValleys) {
    // Higher band(s)
  }

  if (m_nonParabolic) {
    // Non-parabolicity parameter
    double alpha = 0.;
    if (band < m_cb[0].nValleys) {
      // X valley
      alpha = m_cb[0].alpha;
    } else if (band < m_cb[0].nValleys + m_cb[1].nValleys) {
      // L valley
      alpha = m_cb[1].alpha;
    }

    const double p2 = 0.5 * (px * px / mx + py * py / my + pz * pz / mz);
    if (alpha > 0.) {
      const double e = 0.5 * (sqrt(1. + 4 * alpha * p2) - 1.) / alpha;
      const double a = SpeedOfLight / (1. + 2 * alpha * e);
      vx = a * px / mx;
      vy = a * py / my;
      vz = a * pz / mz;
      return e0 + e;
    }
  }

  const double e = 0.5 * (px * px / mx + py * py / my + pz * pz / mz);
  vx = SpeedOfLight * px / mx;
  vy = SpeedOfLight * py / my;
  vz = SpeedOfLight * pz / mz;
  return e0 + e;
}

void MediumSilicon::GetElectronMomentum(const double e, double& px, double& py,
                                        double& pz, int& band) {
  // If the band index is out of range, choose one at random.
  if (band < 0 || band > m_cb[0].nValleys + m_cb[1].nValleys ||
      (e < m_cb[1].eMin && band >= m_cb[0].nValleys) ||
      (e < m_cb[2].eMin && band == m_cb[0].nValleys + m_cb[1].nValleys)) {
    if (e < m_cb[1].eMin) {
      band = int(m_cb[0].nValleys * RndmUniform());
      if (band >= m_cb[0].nValleys) band = m_cb[0].nValleys - 1;
    } else {
      const double dosX = ConductionBandDOS(e, 0);
      const double dosL = ConductionBandDOS(e, m_cb[0].nValleys);
      const double dosG = ConductionBandDOS(e, m_cb[0].nValleys + m_cb[1].nValleys);
      const double dosSum = m_cb[0].nValleys * dosX + m_cb[1].nValleys * dosL + dosG;
      if (dosSum < Small) {
        band = m_cb[0].nValleys + m_cb[1].nValleys;
      } else {
        const double r = RndmUniform() * dosSum;
        if (r < dosX) {
          band = int(m_cb[0].nValleys * RndmUniform());
          if (band >= m_cb[0].nValleys) band = m_cb[0].nValleys - 1;
        } else if (r < dosX + dosL) {
          band = int(m_cb[1].nValleys * RndmUniform());
          if (band >= m_cb[1].nValleys) band = m_cb[1].nValleys - 1;
          band += m_cb[0].nValleys;
        } else {
          band = m_cb[0].nValleys + m_cb[1].nValleys;
        }
      }
    }
    if (m_debug) {
      std::cout << m_className << "::GetElectronMomentum:\n"
                << "    Randomised band index: " << band << "\n";
    }
  }
  const auto k = m_cbIndex[band];
  if (k == 0) {
    // X valleys
    double pstar = sqrt(2. * ElectronMass * e);
    if (m_nonParabolic) {
      pstar *= sqrt(1. + m_cb[0].alpha * e);
    }

    const double ctheta = 1. - 2. * RndmUniform();
    const double stheta = sqrt(1. - ctheta * ctheta);
    const double phi = TwoPi * RndmUniform();

    if (m_anisotropic) {
      const double pl = pstar * sqrt(m_cb[0].mL);
      const double pt = pstar * sqrt(m_cb[0].mT);
      switch (band) {
        case 0:
        case 1:
          // 100
          px = pl * ctheta;
          py = pt * cos(phi) * stheta;
          pz = pt * sin(phi) * stheta;
          break;
        case 2:
        case 3:
          // 010
          px = pt * sin(phi) * stheta;
          py = pl * ctheta;
          pz = pt * cos(phi) * stheta;
          break;
        case 4:
        case 5:
          // 001
          px = pt * cos(phi) * stheta;
          py = pt * sin(phi) * stheta;
          pz = pl * ctheta;
          break;
        default:
          // Other band; should not occur.
          std::cerr << m_className << "::GetElectronMomentum:\n"
                    << "    Unexpected band index (" << band << ").\n";
          px = pstar * stheta * cos(phi);
          py = pstar * stheta * sin(phi);
          pz = pstar * ctheta;
          break;
      }
    } else {
      pstar *= sqrt(m_cb[0].mC);
      px = pstar * cos(phi) * stheta;
      py = pstar * sin(phi) * stheta;
      pz = pstar * ctheta;
    }
  } else if (k == 1) {
    // L valleys
    double pstar = sqrt(2. * ElectronMass * (e - m_cb[1].eMin));
    if (m_nonParabolic) {
      pstar *= sqrt(1. + m_cb[1].alpha * (e - m_cb[1].eMin));
    }
    pstar *= sqrt(m_cb[1].mC);
    RndmDirection(px, py, pz, pstar);
  } else if (k == 2) {
    // Higher band
    const double pstar = sqrt(2. * ElectronMass * e);
    RndmDirection(px, py, pz, pstar);
  }
}

double MediumSilicon::GetElectronNullCollisionRate(const int band) {
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::GetElectronNullCollisionRate:\n"
                << "    Error calculating the collision rates table.\n";
      return 0.;
    }
    m_isChanged = false;
  }
  if (band < 0 || band >= m_cbIndex.size()) {
    std::cerr << m_className << "::GetElectronNullCollisionRate:\n"
              << "    Band index (" << band << ") out of range.\n";
    return 0.;
  }
  return m_cb[m_cbIndex[band]].cfNull;
}

double MediumSilicon::GetElectronCollisionRate(const double e, const int band) {
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n"
              << "    Electron energy must be positive.\n";
    return 0.;
  }

  if (e > m_cb[2].eFinal) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n"
              << "    Collision rate at " << e << " eV (band " << band
              << ") is not included in the current table.\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  }

  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::GetElectronCollisionRate:\n"
                << "    Error calculating the collision rates table.\n";
      return 0.;
    }
    m_isChanged = false;
  }

  if (band < 0 || band >= m_cbIndex.size()) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n"
              << "    Band index (" << band << ") out of range.\n";
    return 0.;
  }
  const size_t k = m_cbIndex[band];
  int iE = int(e / m_cb[k].eStep);
  if (iE >= m_cb[k].nEnergySteps) {
    iE = m_cb[k].nEnergySteps - 1;
  } else if (iE < m_cb[k].iMin) {
    iE = m_cb[k].iMin;
  }
  return m_cb[k].cfTot[iE];
}

bool MediumSilicon::ElectronCollision(const double e, int& type, 
    int& level, double& e1, double& px, double& py, double& pz, 
    std::vector<Secondary>& secondaries, int& band) {
  if (e > m_cb[2].eFinal) {
    std::cerr << m_className << "::ElectronCollision:\n"
              << "    Requested electron energy (" << e << " eV) exceeds the "
              << "current energy range (" << m_cb[2].eFinal << " eV).\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  } else if (e <= 0.) {
    std::cerr << m_className << "::ElectronCollision:\n"
              << "    Electron energy must be greater than zero.\n";
    return false;
  }

  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::ElectronCollision:\n"
                << "    Error calculating the collision rates table.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (band < 0 || band >= m_cbIndex.size()) {
    std::cerr << m_className << "::ElectronCollision:\n"
              << "    Band index (" << band << ") out of range.\n";
    return false;
  }
  const auto k = m_cbIndex[band];
  // Get the energy interval.
  int iE = int(e / m_cb[k].eStep);
  if (iE >= m_cb[k].nEnergySteps) {
    iE = m_cb[k].nEnergySteps - 1;
  } else if (iE < m_cb[k].iMin) {
    iE = m_cb[k].iMin;
  }
  // Select the scattering process.
  const double r = RndmUniform();
  if (r <= m_cb[k].cf[iE][0]) {
    level = 0;
  } else if (r >= m_cb[k].cf[iE][m_cb[k].nLevels - 1]) {
    level = m_cb[k].nLevels - 1;
  } else {
    const auto begin = m_cb[k].cf[iE].cbegin();
    level = std::lower_bound(begin, begin + m_cb[k].nLevels, r) - begin;
  }
  // Get the collision type.
  type = m_cb[k].scatType[level];
  // Fill the collision counters.
  if (k == 0) {
    ++m_nCollElectronDetailed[level];
  } else if (k == 1) {
    ++m_nCollElectronDetailed[m_cb[0].nLevels + level];
  } else if (k == 2) {
    ++m_nCollElectronDetailed[m_cb[0].nLevels + m_cb[1].nLevels + level];
  }
  ++m_nCollElectronBand[band];
  if (type == ElectronCollisionTypeAcousticPhonon) {
    ++m_nCollElectronAcoustic;
  } else if (type == ElectronCollisionTypeOpticalPhonon) {
    ++m_nCollElectronOptical;
  } else if (type == ElectronCollisionTypeIntervalleyG) {
    // Intervalley scattering (g type).
    ++m_nCollElectronIntervalley;
    if (k == 0) {
     // XX. Final valley is in opposite direction.
      switch (band) {
        case 0:
          band = 1;
          break;
        case 1:
          band = 0;
          break;
        case 2:
          band = 3;
          break;
        case 3:
          band = 2;
          break;
        case 4:
          band = 5;
          break;
        case 5:
          band = 4;
          break;
        default:
          break;
      }
    } else if (k == 1) {
      // LL. Randomise the final valley.
      band = int(RndmUniform() * m_cb[1].nValleys);
      if (band >= m_cb[1].nValleys) band = m_cb[1].nValleys;
      band += m_cb[0].nValleys;
    }
  } else if (type == ElectronCollisionTypeIntervalleyF) {
    // Intervalley scattering (f type).
    ++m_nCollElectronIntervalley;
    if (k == 0) {
      // XX. Final valley is perpendicular.
      switch (band) {
        case 0:
        case 1:
          band = int(RndmUniform() * 4) + 2;
          break;
        case 2:
        case 3:
          band = int(RndmUniform() * 4);
          if (band > 1) band += 2;
          break;
        case 4:
        case 5:
          band = int(RndmUniform() * 4);
          break;
      }
    } else if (k == 1) {
      // LL. Randomise the final valley.
      band = int(RndmUniform() * m_cb[1].nValleys);
      if (band >= m_cb[1].nValleys) band = m_cb[1].nValleys;
      band += m_cb[0].nValleys;
    }
  } else if (type == ElectronCollisionTypeInterbandXL) {
    ++m_nCollElectronIntervalley;
    if (k == 0) {
      // XL scattering. Final valley is in L band.
      band = int(RndmUniform() * m_cb[1].nValleys);
      if (band >= m_cb[1].nValleys) band = m_cb[1].nValleys - 1;
      band += m_cb[0].nValleys;
    } else if (k == 1) {
      // LX scattering. Randomise the final valley.
      band = int(RndmUniform() * m_cb[0].nValleys);
      if (band >= m_cb[0].nValleys) band = m_cb[0].nValleys - 1;
    }
  } else if (type == ElectronCollisionTypeInterbandXG) {
    ++m_nCollElectronIntervalley;
    if (k == 0) {
      // XG scattering.
      band = m_cb[0].nValleys + m_cb[1].nValleys;
    } else if (k == 2) {
      // GX scattering. Randomise the final valley.
      band = int(RndmUniform() * m_cb[0].nValleys);
      if (band >= m_cb[0].nValleys) band = m_cb[0].nValleys - 1;
    }
  } else if (type == ElectronCollisionTypeInterbandLG) {
    ++m_nCollElectronIntervalley;
    if (k == 1) {
      // LG scattering
      band = m_cb[0].nValleys + m_cb[1].nValleys;
    } else if (k == 2) {
      // GL scattering. Randomise the final valley.
      band = int(RndmUniform() * m_cb[1].nValleys);
      if (band >= m_cb[1].nValleys) band = m_cb[1].nValleys - 1;
      band += m_cb[0].nValleys;
    }
  } else if (type == ElectronCollisionTypeImpurity) {
    ++m_nCollElectronImpurity;
  } else if (type == ElectronCollisionTypeIonisation) {
    ++m_nCollElectronIonisation;
  } else {
    std::cerr << m_className << "::ElectronCollision:\n"
              << "    Unexpected collision type (" << type << ").\n";
  }

  // Get the energy loss.
  double loss = m_cb[k].energyLoss[level];

  // Ionising collision
  if (type == ElectronCollisionTypeIonisation) {
    double ee = 0., eh = 0.;
    ComputeSecondaries(e, ee, eh);
    loss = ee + eh + m_bandGap;
    // Add the secondary electron.
    Secondary esec;
    esec.type = Particle::Electron;
    esec.energy = ee;
    secondaries.push_back(std::move(esec));
    // Add the hole.
    Secondary hsec;
    hsec.type = Particle::Hole;
    hsec.energy = eh;
    secondaries.push_back(std::move(hsec));
  }

  if (e < loss) loss = e - 0.00001;
  // Update the energy.
  e1 = e - loss;
  if (e1 < Small) e1 = Small;

  // Update the momentum.
  if (k == 0) {
    // X valleys
    double pstar = sqrt(2. * ElectronMass * e1);
    if (m_nonParabolic) {
      pstar *= sqrt(1. + m_cb[0].alpha * e1);
    }

    const double ctheta = 1. - 2. * RndmUniform();
    const double stheta = sqrt(1. - ctheta * ctheta);
    const double phi = TwoPi * RndmUniform();

    if (m_anisotropic) {
      const double pl = pstar * sqrt(m_cb[0].mL);
      const double pt = pstar * sqrt(m_cb[0].mT);
      switch (band) {
        case 0:
        case 1:
          // 100
          px = pl * ctheta;
          py = pt * cos(phi) * stheta;
          pz = pt * sin(phi) * stheta;
          break;
        case 2:
        case 3:
          // 010
          px = pt * sin(phi) * stheta;
          py = pl * ctheta;
          pz = pt * cos(phi) * stheta;
          break;
        case 4:
        case 5:
          // 001
          px = pt * cos(phi) * stheta;
          py = pt * sin(phi) * stheta;
          pz = pl * ctheta;
          break;
        default:
          return false;
      }
    } else {
      pstar *= sqrt(m_cb[0].mC);
      px = pstar * cos(phi) * stheta;
      py = pstar * sin(phi) * stheta;
      pz = pstar * ctheta;
    }
    return true;

  } else if (k == 1) {
    // L valleys
    double pstar = sqrt(2. * ElectronMass * (e1 - m_cb[1].eMin));
    if (m_nonParabolic) {
      pstar *= sqrt(1. + m_cb[1].alpha * (e1 - m_cb[1].eMin));
    }
    pstar *= sqrt(m_cb[1].mC);
    RndmDirection(px, py, pz, pstar);
    return true;
  }
  const double pstar = sqrt(2. * ElectronMass * e1);
  RndmDirection(px, py, pz, pstar);
  return true;
}

void MediumSilicon::ResetCollisionCounters() {
  m_nCollElectronAcoustic = m_nCollElectronOptical = 0;
  m_nCollElectronIntervalley = 0;
  m_nCollElectronImpurity = 0;
  m_nCollElectronIonisation = 0;
  const auto nLevels = m_cb[0].nLevels + m_cb[1].nLevels + m_cb[2].nLevels;
  m_nCollElectronDetailed.assign(nLevels, 0);
  const auto nBands = m_cb[0].nValleys + m_cb[1].nValleys + 1;
  m_nCollElectronBand.assign(nBands, 0);
}

unsigned int MediumSilicon::GetNumberOfElectronCollisions() const {
  return m_nCollElectronAcoustic + m_nCollElectronOptical +
         m_nCollElectronIntervalley + m_nCollElectronImpurity +
         m_nCollElectronIonisation;
}

unsigned int MediumSilicon::GetNumberOfLevels() const {
  return m_cb[0].nLevels + m_cb[1].nLevels + m_cb[2].nLevels;
}

unsigned int MediumSilicon::GetNumberOfElectronCollisions(
    const unsigned int level) const {
  const auto nLevels = m_cb[0].nLevels + m_cb[1].nLevels + m_cb[2].nLevels;
  if (level >= nLevels) {
    std::cerr << m_className << "::GetNumberOfElectronCollisions:\n"
              << "    Scattering rate term (" << level << ") does not exist.\n";
    return 0;
  }
  return m_nCollElectronDetailed[level];
}

unsigned int MediumSilicon::GetNumberOfElectronBands() const {
  return m_cb[0].nValleys + m_cb[1].nValleys + 1;
}

int MediumSilicon::GetElectronBandPopulation(const int band) {
  const int nBands = m_cb[0].nValleys + m_cb[1].nValleys + 1;
  if (band < 0 || band >= nBands) {
    std::cerr << m_className << "::GetElectronBandPopulation:\n";
    std::cerr << "    Band index (" << band << ") out of range.\n";
    return 0;
  }
  return m_nCollElectronBand[band];
}

bool MediumSilicon::GetOpticalDataRange(double& emin, double& emax,
                                        const unsigned int i) {
  if (i != 0) {
    std::cerr << m_className << "::GetOpticalDataRange: Index out of range.\n";
  }

  // Make sure the optical data table has been loaded.
  if (m_opticalDataEnergies.empty()) {
    if (!LoadOpticalData(m_opticalDataFile)) {
      std::cerr << m_className << "::GetOpticalDataRange:\n"
                << "    Optical data table could not be loaded.\n";
      return false;
    }
  }

  emin = m_opticalDataEnergies.front();
  emax = m_opticalDataEnergies.back();
  if (m_debug) {
    std::cout << m_className << "::GetOpticalDataRange:\n"
              << "    " << emin << " < E [eV] < " << emax << "\n";
  }
  return true;
}

bool MediumSilicon::GetDielectricFunction(const double e, double& eps1,
                                          double& eps2, const unsigned int i) {
  if (i != 0) {
    std::cerr << m_className + "::GetDielectricFunction: Index out of range.\n";
    return false;
  }

  // Make sure the optical data table has been loaded.
  if (m_opticalDataEnergies.empty()) {
    if (!LoadOpticalData(m_opticalDataFile)) {
      std::cerr << m_className << "::GetDielectricFunction:\n";
      std::cerr << "    Optical data table could not be loaded.\n";
      return false;
    }
  }

  // Make sure the requested energy is within the range of the table.
  const double emin = m_opticalDataEnergies.front();
  const double emax = m_opticalDataEnergies.back();
  if (e < emin || e > emax) {
    std::cerr << m_className << "::GetDielectricFunction:\n"
              << "    Requested energy (" << e << " eV) "
              << " is outside the range of the optical data table.\n"
              << "    " << emin << " < E [eV] < " << emax << "\n";
    eps1 = eps2 = 0.;
    return false;
  }

  // Locate the requested energy in the table.
  const auto begin = m_opticalDataEnergies.cbegin();
  const auto it1 = std::upper_bound(begin, m_opticalDataEnergies.cend(), e);
  if (it1 == begin) {
    eps1 = m_opticalDataEpsilon.front().first;
    eps2 = m_opticalDataEpsilon.front().second;
    return true;
  }
  const auto it0 = std::prev(it1);

  // Interpolate the real part of dielectric function.
  const double x0 = *it0;
  const double x1 = *it1;
  const double lnx0 = log(*it0);
  const double lnx1 = log(*it1);
  const double lnx = log(e);
  const double y0 = m_opticalDataEpsilon[it0 - begin].first;
  const double y1 = m_opticalDataEpsilon[it1 - begin].first;
  if (y0 <= 0. || y1 <= 0.) {
    // Use linear interpolation if one of the values is negative.
    eps1 = y0 + (e - x0) * (y1 - y0) / (x1 - x0);
  } else {
    // Otherwise use log-log interpolation.
    const double lny0 = log(y0);
    const double lny1 = log(y1);
    eps1 = lny0 + (lnx - lnx0) * (lny1 - lny0) / (lnx1 - lnx0);
    eps1 = exp(eps1);
  }

  // Interpolate the imaginary part of dielectric function,
  // using log-log interpolation.
  const double lnz0 = log(m_opticalDataEpsilon[it0 - begin].second);
  const double lnz1 = log(m_opticalDataEpsilon[it1 - begin].second);
  eps2 = lnz0 + (lnx - lnx0) * (lnz1 - lnz0) / (lnx1 - lnx0);
  eps2 = exp(eps2);
  return true;
}

bool MediumSilicon::Initialise() {
  if (!m_isChanged) {
    if (m_debug) {
      std::cerr << m_className << "::Initialise: Nothing changed.\n";
    }
    return true;
  }
  if (!UpdateTransportParameters()) {
    std::cerr << m_className << "::Initialise:\n    Error preparing "
              << "transport parameters/calculating collision rates.\n";
    return false;
  }
  return true;
}

bool MediumSilicon::UpdateTransportParameters() {
  std::lock_guard<std::mutex> guard(m_mutex);

  // Calculate impact ionisation coefficients.
  UpdateImpactIonisation();

  if (!m_hasUserMobility) {
    // Calculate lattice mobility.
    UpdateLatticeMobility();

    // Calculate doping mobility
    switch (m_dopingMobilityModel) {
      case DopingMobility::Minimos:
        UpdateDopingMobilityMinimos();
        break;
      case DopingMobility::Masetti:
        UpdateDopingMobilityMasetti();
        break;
      default:
        std::cerr << m_className << "::UpdateTransportParameters:\n    "
                  << "Unknown doping mobility model. Program bug!\n";
        break;
    }
  }

  // Calculate saturation velocity
  if (!m_hasUserSaturationVelocity) {
    UpdateSaturationVelocity();
  }

  // Calculate high field saturation parameters
  if (m_highFieldMobilityModel == HighFieldMobility::Canali) {
    UpdateHighFieldMobilityCanali();
  }

  if (m_debug) {
    std::cout << m_className << "::UpdateTransportParameters:\n"
              << "    Low-field mobility [cm2 V-1 ns-1]\n"
              << "      Electrons: " << m_eMobility << "\n"
              << "      Holes:     " << m_hMobility << "\n";
    if (m_highFieldMobilityModel == HighFieldMobility::Constant) {
      std::cout << "    Mobility is not field-dependent.\n";
    } else {
      std::cout << "    Saturation velocity [cm / ns]\n"
                << "      Electrons: " << m_eSatVel << "\n"
                << "      Holes:     " << m_hSatVel << "\n";
    }
  }

  if (!ElectronScatteringRates()) return false;
  if (!HoleScatteringRates()) return false;

  // Reset the collision counters.
  ResetCollisionCounters();

  return true;
}

void MediumSilicon::UpdateLatticeMobility() {

  // Temperature normalized to 300 K
  const double t = m_temperature / 300.;

  switch (m_latticeMobilityModel) {
    case LatticeMobility::Minimos:
      // - S. Selberherr, W. Haensch, M. Seavey, J. Slotboom,
      //   Solid State Electronics 33 (1990), 1425
      // - Minimos 6.1 User's Guide (1999)
      m_eLatticeMobility = 1.43e-6 * pow(t, -2.);
      m_hLatticeMobility = 0.46e-6 * pow(t, -2.18);
      break;
    case LatticeMobility::Sentaurus:
      // - C. Lombardi et al., IEEE Trans. CAD 7 (1988), 1164
      // - Sentaurus Device User Guide (2007)
      m_eLatticeMobility = 1.417e-6 * pow(t, -2.5);
      m_hLatticeMobility = 0.4705e-6 * pow(t, -2.2);
      break;
    case LatticeMobility::Reggiani:
      // M. A. Omar, L. Reggiani, Solid State Electronics 30 (1987), 693
      m_eLatticeMobility = 1.320e-6 * pow(t, -2.);
      m_hLatticeMobility = 0.460e-6 * pow(t, -2.2);
      break;
    default:
      std::cerr << m_className << "::UpdateLatticeMobility:\n"
                << "    Unknown lattice mobility model. Program bug!\n";
      break;
  }
}

void MediumSilicon::UpdateDopingMobilityMinimos() {
  // References:
  // - S. Selberherr, W. Haensch, M. Seavey, J. Slotboom,
  //   Solid State Electronics 33 (1990), 1425-1436
  // - Minimos 6.1 User's Guide (1999)

  // Mobility reduction due to ionised impurity scattering
  // Surface term not taken into account
  double eMuMin = 0.080e-6;
  double hMuMin = 0.045e-6;
  if (m_temperature > 200.) {
    const double c0 = pow(m_temperature / 300., -0.45);
    eMuMin *= c0;
    hMuMin *= c0;
  } else {
    const double c0 = pow(2. / 3., -0.45) * pow(m_temperature / 200., -0.15);
    eMuMin *= c0;
    hMuMin *= c0;
  }
  const double c1 = pow(m_temperature / 300., 3.2);
  const double eRefC = 1.12e17 * c1;
  const double hRefC = 2.23e17 * c1;
  const double alpha = 0.72 * pow(m_temperature / 300., 0.065);
  // Assume impurity concentration equal to doping concentration
  m_eMobility = eMuMin +
                (m_eLatticeMobility - eMuMin) /
                    (1. + pow(m_dopingConcentration / eRefC, alpha));
  m_hMobility = hMuMin +
                (m_hLatticeMobility - hMuMin) /
                    (1. + pow(m_dopingConcentration / hRefC, alpha));
}

void MediumSilicon::UpdateDopingMobilityMasetti() {
  // Reference:
  // - G. Masetti, M. Severi, S. Solmi,
  //   IEEE Trans. Electron Devices 30 (1983), 764-769
  // - Sentaurus Device User Guide (2007)
  // - Minimos NT User Guide (2004)

  if (m_dopingConcentration < 1.e13) {
    m_eMobility = m_eLatticeMobility;
    m_hMobility = m_hLatticeMobility;
    return;
  }

  // Parameters adopted from Minimos NT User Guide
  constexpr double eMuMin1 = 0.0522e-6;
  constexpr double eMuMin2 = 0.0522e-6;
  constexpr double eMu1 = 0.0434e-6;
  constexpr double hMuMin1 = 0.0449e-6;
  constexpr double hMuMin2 = 0.;
  constexpr double hMu1 = 0.029e-6;
  constexpr double eCr = 9.68e16;
  constexpr double eCs = 3.42e20;
  constexpr double hCr = 2.23e17;
  constexpr double hCs = 6.10e20;
  constexpr double hPc = 9.23e16;
  constexpr double eAlpha = 0.68;
  constexpr double eBeta = 2.;
  constexpr double hAlpha = 0.719;
  constexpr double hBeta = 2.;

  m_eMobility = eMuMin1 +
                (m_eLatticeMobility - eMuMin2) /
                    (1. + pow(m_dopingConcentration / eCr, eAlpha)) -
                eMu1 / (1. + pow(eCs / m_dopingConcentration, eBeta));

  m_hMobility = hMuMin1 * exp(-hPc / m_dopingConcentration) +
                (m_hLatticeMobility - hMuMin2) /
                    (1. + pow(m_dopingConcentration / hCr, hAlpha)) -
                hMu1 / (1. + pow(hCs / m_dopingConcentration, hBeta));
}

void MediumSilicon::UpdateSaturationVelocity() {

  switch (m_saturationVelocityModel) {
    case SaturationVelocity::Minimos:
      // - R. Quay, C. Moglestue, V. Palankovski, S. Selberherr,
      //   Materials Science in Semiconductor Processing 3 (2000), 149
      // - Minimos NT User Guide (2004)
      m_eSatVel = 1.e-2 / (1. + 0.74 * (m_temperature / 300. - 1.));
      m_hSatVel = 0.704e-2 / (1. + 0.37 * (m_temperature / 300. - 1.));
      break;
    case SaturationVelocity::Canali:
      // - C. Canali, G. Majni, R. Minder, G. Ottaviani,
      //   IEEE Transactions on Electron Devices 22 (1975), 1045
      // - Sentaurus Device User Guide (2007)
      m_eSatVel = 1.07e-2 * pow(300. / m_temperature, 0.87);
      m_hSatVel = 8.37e-3 * pow(300. / m_temperature, 0.52);
      break;
    case SaturationVelocity::Reggiani:
      // M. A. Omar, L. Reggiani, Solid State Electronics 30 (1987), 693
      m_eSatVel = 1.470e-2 * sqrt(tanh(150. / m_temperature));
      m_hSatVel = 0.916e-2 * sqrt(tanh(300. / m_temperature));
      break;
    default:
      std::cerr << m_className << "::UpdateSaturationVelocity:\n" 
                << "    Unknown saturation velocity model. Program bug!\n";
      break;
  }
}

void MediumSilicon::UpdateHighFieldMobilityCanali() {
  // References:
  // - C. Canali, G. Majni, R. Minder, G. Ottaviani,
  //   IEEE Transactions on Electron Devices 22 (1975), 1045-1047
  // - Sentaurus Device User Guide (2007)

  // Temperature dependent exponent in high-field mobility formula
  m_eBetaCanali = 1.109 * pow(m_temperature / 300., 0.66);
  m_hBetaCanali = 1.213 * pow(m_temperature / 300., 0.17);
  m_eBetaCanaliInv = 1. / m_eBetaCanali;
  m_hBetaCanaliInv = 1. / m_hBetaCanali;
}

void MediumSilicon::UpdateImpactIonisation() {

  if (m_impactIonisationModel == ImpactIonisation::VanOverstraeten ||
      m_impactIonisationModel == ImpactIonisation::Grant) {

    // Temperature dependence as in Sentaurus Device
    // Optical phonon energy
    constexpr double hbarOmega = 0.063;
    // Temperature scaling coefficient
    const double gamma =
        tanh(hbarOmega / (2. * BoltzmannConstant * 300.)) /
        tanh(hbarOmega / (2. * BoltzmannConstant * m_temperature));

    if (m_impactIonisationModel == ImpactIonisation::VanOverstraeten) {
      // - R. van Overstraeten and H. de Man,
      //   Solid State Electronics 13 (1970), 583
      // - W. Maes, K. de Meyer and R. van Overstraeten,
      //   Solid State Electronics 33 (1990), 705
      // - Sentaurus Device User Guide (2016)

      // Low field coefficients taken from Maes, de Meyer, van Overstraeten
      // eImpactA0 = gamma * 3.318e5;
      // eImpactB0 = gamma * 1.135e6;
      m_eImpactA0 = gamma * 7.03e5;
      m_eImpactB0 = gamma * 1.231e6;
      m_eImpactA1 = gamma * 7.03e5;
      m_eImpactB1 = gamma * 1.231e6;

      m_hImpactA0 = gamma * 1.582e6;
      m_hImpactB0 = gamma * 2.036e6;
      m_hImpactA1 = gamma * 6.71e5;
      m_hImpactB1 = gamma * 1.693e6;
    } else {
      // W. N. Grant, Solid State Electronics 16 (1973), 1189
      // Sentaurus Device User Guide (2007)
      m_eImpactA0 = 2.60e6 * gamma;
      m_eImpactB0 = 1.43e6 * gamma;
      m_eImpactA1 = 6.20e5 * gamma;
      m_eImpactB1 = 1.08e6 * gamma;
      m_eImpactA2 = 5.05e5 * gamma;
      m_eImpactB2 = 9.90e5 * gamma;

      m_hImpactA0 = 2.00e6 * gamma;
      m_hImpactB0 = 1.97e6 * gamma;
      m_hImpactA1 = 5.60e5 * gamma;
      m_hImpactB1 = 1.32e6 * gamma;
    }
  } else if (m_impactIonisationModel == ImpactIonisation::Massey) {
    // D. J. Massey, J. P. R. David, and G. J. Rees,
    // IEEE Trans. Electron Devices 53 (2006), 2328
    m_eImpactA0 = 4.43e5;
    m_eImpactB0 = 9.66e5 + 4.99e2 * m_temperature;

    m_hImpactA0 = 1.13e6;
    m_hImpactB0 = 1.71e6 + 1.09e3 * m_temperature;
  } else if (m_impactIonisationModel == ImpactIonisation::Okuto) {
    const double dt = m_temperature - 300.;
    m_eImpactA0 = 0.426 * (1. + 3.05e-4 * dt);
    m_hImpactA0 = 0.243 * (1. + 5.35e-4 * dt);
    m_eImpactB0 = 4.81e5 * (1. + 6.86e-4 * dt);
    m_hImpactB0 = 6.53e5 * (1. + 5.67e-4 * dt);
  } else {
    std::cerr << m_className << "::UpdateImpactIonisation:\n"
              << "    Unknown impact ionisation model. Program bug!\n";
  }
}

double MediumSilicon::ElectronMobility(const double emag) const {

  if (emag < Small) return 0.;

  if (m_highFieldMobilityModel == HighFieldMobility::Minimos) {
    // Minimos User's Guide (1999)
    const double r = 2 * m_eMobility * emag / m_eSatVel;
    return 2. * m_eMobility / (1. + sqrt(1. + r * r));
  } else if (m_highFieldMobilityModel == HighFieldMobility::Canali) {
    // Sentaurus Device User Guide (2007)
    const double r = m_eMobility * emag / m_eSatVel;
    return m_eMobility / pow(1. + pow(r, m_eBetaCanali), m_eBetaCanaliInv);
  } else if (m_highFieldMobilityModel == HighFieldMobility::Reggiani) {
    // M. A. Omar, L. Reggiani, Solid State Electronics 30 (1987), 693
    const double r = m_eMobility * emag / m_eSatVel;
    constexpr double k = 1. / 1.5;
    return m_eMobility / pow(1. + pow(r, 1.5), k);
  }
  return m_eMobility;
}

double MediumSilicon::ElectronAlpha(const double emag) const {

  if (emag < Small) return 0.;
  
  if (m_impactIonisationModel == ImpactIonisation::VanOverstraeten) {
    // - R. van Overstraeten and H. de Man,
    //   Solid State Electronics 13 (1970), 583
    // - W. Maes, K. de Meyer and R. van Overstraeten,
    //   Solid State Electronics 33 (1990), 705
    // - Sentaurus Device User Guide (2016)
    if (emag < 4e5) {
      return m_eImpactA0 * exp(-m_eImpactB0 / emag);
    } else {
      return m_eImpactA1 * exp(-m_eImpactB1 / emag);
    }
  } else if (m_impactIonisationModel == ImpactIonisation::Grant) {
    // W. N. Grant, Solid State Electronics 16 (1973), 1189
    if (emag < 2.4e5) {
      return m_eImpactA0 * exp(-m_eImpactB0 / emag);
    } else if (emag < 5.3e5) {
      return m_eImpactA1 * exp(-m_eImpactB1 / emag);
    } else {
      return m_eImpactA2 * exp(-m_eImpactB2 / emag);
    }
  } else if (m_impactIonisationModel == ImpactIonisation::Massey) {
    return m_eImpactA0 * exp(-m_eImpactB0 / emag);
  } else if (m_impactIonisationModel == ImpactIonisation::Okuto) {
    const double f = m_eImpactB0 / emag;
    return m_eImpactA0 * emag * exp(-f * f);
  }
  std::cerr << m_className << "::ElectronAlpha: Unknown model. Program bug!\n";
  return 0.;
}

double MediumSilicon::HoleMobility(const double emag) const {

  if (emag < Small) return 0.;

  if (m_highFieldMobilityModel == HighFieldMobility::Minimos) {
    // Minimos User's Guide (1999)
    return m_hMobility / (1. + m_hMobility * emag / m_eSatVel);
  } else if (m_highFieldMobilityModel == HighFieldMobility::Canali) {
    // Sentaurus Device User Guide (2007)
    const double r = m_hMobility * emag / m_hSatVel;
    return m_hMobility / pow(1. + pow(r, m_hBetaCanali), m_hBetaCanaliInv);
  } else if (m_highFieldMobilityModel == HighFieldMobility::Reggiani) {
    // M. A. Omar, L. Reggiani, Solid State Electronics 30 (1987), 693
    const double r = m_hMobility * emag / m_hSatVel;
    return m_hMobility / sqrt(1. + r * r);
  }
  return m_hMobility;
}

double MediumSilicon::HoleAlpha(const double emag) const {

  if (emag < Small) return 0.;

  if (m_impactIonisationModel == ImpactIonisation::VanOverstraeten) {
    // - R. van Overstraeten and H. de Man,
    //   Solid State Electronics 13 (1970), 583
    // - Sentaurus Device User Guide (2016)
    if (emag < 4e5) {
      return m_hImpactA0 * exp(-m_hImpactB0 / emag);
    } else {
      return m_hImpactA1 * exp(-m_hImpactB1 / emag);
    }
  } else if (m_impactIonisationModel == ImpactIonisation::Grant) {
    // W. N. Grant, Solid State Electronics 16 (1973), 1189
    if (emag < 5.3e5) {
      return m_hImpactA0 * exp(-m_hImpactB0 / emag);
    } else {
      return m_hImpactA1 * exp(-m_hImpactB1 / emag);
    } 
  } else if (m_impactIonisationModel == ImpactIonisation::Massey) {
    return m_hImpactA0 * exp(-m_hImpactB0 / emag);
  } else if (m_impactIonisationModel == ImpactIonisation::Okuto) {
    const double f = m_hImpactB0 / emag;
    return m_hImpactA0 * emag * exp(-f * f);
  }
  std::cerr << m_className << "::HoleAlpha: Unknown model. Program bug!\n";
  return 0.;
}

bool MediumSilicon::LoadOpticalData(const std::string& filename) {
  // Clear the optical data table.
  m_opticalDataEnergies.clear();
  m_opticalDataEpsilon.clear();

  // Get the path to the data directory.
  char* pPath = getenv("GARFIELD_HOME");
  if (pPath == 0) {
    std::cerr << m_className << "::LoadOpticalData:\n"
              << "    Environment variable GARFIELD_HOME is not set.\n";
    return false;
  }
  const std::string filepath = std::string(pPath) + "/Data/" + filename;

  // Open the file.
  std::ifstream infile(filepath);
  // Make sure the file could actually be opened.
  if (!infile) {
    std::cerr << m_className << "::LoadOpticalData:\n"
              << "    Error opening file " << filename << ".\n";
    return false;
  }

  double lastEnergy = -1.;
  double energy, eps1, eps2, loss;
  // Read the file line by line.
  std::string line;
  std::istringstream dataStream;
  int i = 0;
  while (!infile.eof()) {
    ++i;
    // Read the next line.
    std::getline(infile, line);
    // Strip white space from the beginning of the line.
    ltrim(line);
    if (line.empty()) continue;
    // Skip comments.
    if (line[0] == '#' || line[0] == '*' || (line[0] == '/' && line[1] == '/'))
      continue;
    // Extract the values.
    dataStream.str(line);
    dataStream >> energy >> eps1 >> eps2 >> loss;
    if (dataStream.eof()) break;
    // Check if the data has been read correctly.
    if (infile.fail()) {
      std::cerr << m_className << "::LoadOpticalData:\n    Error reading file "
                << filename << " (line " << i << ").\n";
      return false;
    }
    // Reset the stringstream.
    dataStream.str("");
    dataStream.clear();
    // Make sure the values make sense.
    // The table has to be in ascending order
    //  with respect to the photon energy.
    if (energy <= lastEnergy) {
      std::cerr << m_className << "::LoadOpticalData:\n    Table is not in "
                << "monotonically increasing order (line " << i << ").\n"
                << "    " << lastEnergy << "  " << energy << "  " << eps1
                << "  " << eps2 << "\n";
      return false;
    }
    // The imaginary part of the dielectric function has to be positive.
    if (eps2 < 0.) {
      std::cerr << m_className << "::LoadOpticalData:\n    Negative value "
                << "of the loss function (line " << i << ").\n";
      return false;
    }
    // Ignore negative photon energies.
    if (energy <= 0.) continue;
    // Add the values to the list.
    m_opticalDataEnergies.emplace_back(energy);
    m_opticalDataEpsilon.emplace_back(std::make_pair(eps1, eps2));
    lastEnergy = energy;
  }

  if (m_opticalDataEnergies.empty()) {
    std::cerr << m_className << "::LoadOpticalData:\n"
              << "    Import of data from file " << filepath << "failed.\n"
              << "    No valid data found.\n";
    return false;
  }

  if (m_debug) {
    std::cout << m_className << "::LoadOpticalData:\n    Read "
              << m_opticalDataEnergies.size() << " values from file "
              << filepath << "\n";
  }
  return true;
}

bool MediumSilicon::ElectronScatteringRates() {
  // Reset the scattering rates.
  const size_t nC = m_cb.size(); 
  for (size_t i = 0; i < nC; ++i) {
    m_cb[i].cfTot.assign(m_cb[i].nEnergySteps, 0.);
    m_cb[i].cf.assign(m_cb[i].nEnergySteps, std::vector<double>());
    m_cb[i].energyLoss.clear();
    m_cb[i].scatType.clear();
    m_cb[i].cfNull = 0.;
    m_cb[i].nLevels = 0;
  }
  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  // Fill the scattering rate tables.
  const int kX = 0;
  const int kL = m_cb[0].nValleys;
  const int kG = m_cb[0].nValleys + m_cb[1].nValleys;

  // Acoustic phonon intraband scattering
  // Acoustic deformation potential [eV]
  constexpr double dp = 9.;
  AcousticScatteringRates(rho, kbt, dp, m_cb[0], kX);
  AcousticScatteringRates(rho, kbt, dp, m_cb[1], kL);
  AcousticScatteringRates(rho, kbt, dp, m_cb[2], kG);

  // Coupling constant [eV/cm]
  constexpr double dtk = 2.2e8;
  // Phonon energy [eV]
  constexpr double eph = 63.0e-3;
  // OpticalScatteringRates(rho, kbt, dtk, eph, m_cb[1], kL);
  OpticalScatteringRates(rho, kbt, dtk, eph, m_cb[2], kG);

  ImpurityScatteringRates(kbt, m_cb[0]);
  ImpurityScatteringRates(kbt, m_cb[1]);

  // Intervalley scattering
  // Number of equivalent valleys
  constexpr int zX = 6;
  constexpr int zL = 8;
  constexpr int zG = 1;
  // XX
  // g-type scattering: transition between opposite axes (multiplicity 1)
  // TA (g) - LA (g) - LO (g)
  std::array<double, 3> dXXg = {0.5e8, 0.8e8, 1.1e9}; 
  std::array<double, 3> eXXg = {12.06e-3, 18.53e-3, 62.04e-3};
  for (size_t j = 0; j < 3; ++j) {
    IntervalleyScatteringRates(rho, kbt, dXXg[j], eXXg[j], kX, 1, m_cb[0].eMin,
                               ElectronCollisionTypeIntervalleyG, m_cb[0]);
  }
  // f-type scattering: transition between orthogonal axes (multiplicity 4)
  // TA (f) - LA (f) - TO (f)
  std::array<double, 3> dXXf = {0.3e8, 2.0e8, 2.0e8};
  std::array<double, 3> eXXf = {12.06e-3, 18.53e-3, 62.04e-3};
  for (size_t j = 0; j < 3; ++j) {
    IntervalleyScatteringRates(rho, kbt, dXXf[j], eXXf[j], kX, 4, m_cb[0].eMin,
                               ElectronCollisionTypeIntervalleyF, m_cb[0]);
  }
  // XL
  // - M. Lundstrom, Fundamentals of carrier transport
  // - M. Martin et al.,
  //   Semicond. Sci. Technol. 8, 1291-1297
  std::array<double, 4> dXL = {2.e8, 2.e8, 2.e8, 2.e8};
  std::array<double, 4> eXL = {58.e-3, 55.e-3, 41.e-3, 17.e-3};
  for (size_t j = 0; j < 4; ++j) {
    IntervalleyScatteringRates(rho, kbt, dXL[j], eXL[j], kL, zL, m_cb[1].eMin,
                               ElectronCollisionTypeInterbandXL, m_cb[0]);
    IntervalleyScatteringRates(rho, kbt, dXL[j], eXL[j], kX, zX, m_cb[0].eMin,
                               ElectronCollisionTypeInterbandXL, m_cb[1]);
  }
  // LL
  //  - K. Hess (editor),
  //    Monte Carlo device simulation: full band and beyond
  //    Chapter 5
  //  - M. J. Martin et al.,
  //    Semicond. Sci. Technol. 8, 1291-1297
  constexpr double dLL = 2.63e8;
  constexpr double eLL = 38.87e-3;
  IntervalleyScatteringRates(rho, kbt, dLL, eLL, kL, zL - 1, m_cb[1].eMin, 
                             ElectronCollisionTypeIntervalleyF, m_cb[1]);
  // XG, LG
  // Average of XG and LG
  constexpr double dG = 2.43e8;
  constexpr double eG = 37.65e-3;
  IntervalleyScatteringRates(rho, kbt, dG, eG, kG, zG, m_cb[2].eMin,
                             ElectronCollisionTypeInterbandXG, m_cb[0]);
  IntervalleyScatteringRates(rho, kbt, dG, eG, kG, zG, m_cb[2].eMin, 
                             ElectronCollisionTypeInterbandLG, m_cb[1]);
  IntervalleyScatteringRates(rho, kbt, dG, eG, kX, zX, m_cb[0].eMin, 
                             ElectronCollisionTypeInterbandXG, m_cb[2]);
  IntervalleyScatteringRates(rho, kbt, dG, eG, kL, zL, m_cb[1].eMin, 
                             ElectronCollisionTypeInterbandLG, m_cb[2]);

  // Impact ionisation.

  // - E. Cartier, M. V. Fischetti, E. A. Eklund and F. R. McFeely,
  //   Appl. Phys. Lett 62, 3339-3341
  // - DAMOCLES web page: www.research.ibm.com/DAMOCLES
  // Coefficients [ns-1]
  const std::vector<double> pXL = {6.25e1, 3.e3, 6.8e5};
  // Threshold energies [eV]
  const std::vector<double> ethXL = {1.2, 1.8, 3.45};
  IonisationRates(pXL, ethXL, m_cb[0]);
  IonisationRates(pXL, ethXL, m_cb[1]);

  // - E. Cartier, M. V. Fischetti, E. A. Eklund and F. R. McFeely,
  //   Appl. Phys. Lett 62, 3339-3341
  // - S. Tanuma, C. J. Powell and D. R. Penn
  //   Surf. Interface Anal. (2010)
  // Coefficients [ns-1]
  const std::vector<double> pG = {6.25e1, 3.e3, 6.8e5};
  // Threshold energies [eV]
  const std::vector<double> ethG = {1.2, 1.8, 3.45};
  IonisationRates(pG, ethG, m_cb[2]);

  if (m_debug) {
    std::cout << m_className << "::ElectronScatteringRates:\n"
              << "    " << m_cb[0].nLevels << " X-valley scattering terms\n"
              << "    " << m_cb[1].nLevels << " L-valley scattering terms\n"
              << "    " << m_cb[2].nLevels << " higher band scattering terms\n";
  }

  for (size_t k = 0; k < nC; ++k) {
    std::ofstream outfile;
    if (m_cfOutput) {
      if (k == 0) {
        outfile.open("ratesX.txt", std::ios::out);
      } else if (k == 1) {
        outfile.open("ratesL.txt", std::ios::out);
      } else {
        outfile.open("ratesG.txt", std::ios::out);
      }
    }
    for (int i = 0; i < m_cb[k].nEnergySteps; ++i) {
      // Sum up the scattering rates of all processes.
      for (int j = 0; j < m_cb[k].nLevels; ++j) {
        m_cb[k].cfTot[i] += m_cb[k].cf[i][j];
      }
      if (m_cfOutput) {
        outfile << i * m_cb[k].eStep << " " << m_cb[k].cfTot[i] << " ";
        for (int j = 0; j < m_cb[k].nLevels; ++j) {
          outfile << m_cb[k].cf[i][j] << " ";
        }
        outfile << "\n";
      }
      if (m_cb[k].cfTot[i] > m_cb[k].cfNull) {
        m_cb[k].cfNull = m_cb[k].cfTot[i];
      }
      // Make sure the total scattering rate is positive.
      if (m_cb[k].cfTot[i] <= 0.) {
        std::cerr << m_className << "::ElectronScatteringRates:\n"
                  << "    Scattering rate at " << i * m_cb[k].eStep 
                  << " eV <= 0.\n";
        // outfile.close();
        return false;
      }
      // Normalise the rates.
      for (int j = 0; j < m_cb[k].nLevels; ++j) {
        m_cb[k].cf[i][j] /= m_cb[k].cfTot[i];
        if (j > 0) m_cb[k].cf[i][j] += m_cb[k].cf[i][j - 1];
      }
    }
    if (m_cfOutput) outfile.close();
  }
  return true;
}

bool MediumSilicon::AcousticScatteringRates(
    const double rho, const double kbt, const double dp,
    Band& band, const int k) {

  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705

  // Longitudinal velocity of sound [cm/ns]
  constexpr double u = 9.04e-4;

  // Prefactor for acoustic deformation potential scattering
  const double c = TwoPi * SpeedOfLight * SpeedOfLight * kbt * dp * dp /
                   (Hbar * u * u * rho);

  for (int i = 0; i < band.nEnergySteps; ++i) {
    const double en = Small + i * band.eStep;
    band.cf[i].push_back(c * ConductionBandDOS(en, k));
  }
  
  // Assume that energy loss is negligible.
  band.energyLoss.push_back(0.);
  band.scatType.push_back(ElectronCollisionTypeAcousticPhonon);
  band.nLevels += 1;
  return true;
}

bool MediumSilicon::OpticalScatteringRates(
    const double rho, const double kbt, const double dtk, const double eph,
    Band& band, const int k) {

  // Reference:
  //  - K. Hess (editor),
  //    Monte Carlo device simulation: full band and beyond
  //    Chapter 5
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705
  //  - M. Lundstrom,
  //    Fundamentals of carrier transport

  // Phonon occupation numbers
  const double nocc = 1. / (exp(eph / kbt) - 1);
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c = c0 * dtk * dtk / eph;

  for (int i = 0; i < band.nEnergySteps; ++i) {
    const double en = i * band.eStep;
    // Absorption
    if (en > band.eMin) {
      band.cf[i].push_back(c * nocc * ConductionBandDOS(en + eph, k));
    } else {
      band.cf[i].push_back(0.);
    }
    // Emission
    if (en > band.eMin + eph) {
      band.cf[i].push_back(c * (nocc + 1) * ConductionBandDOS(en - eph, k));
    } else {
      band.cf[i].push_back(0.);
    }
  }

  // Absorption
  band.energyLoss.push_back(-eph);
  // Emission
  band.energyLoss.push_back(eph);
  band.scatType.push_back(ElectronCollisionTypeOpticalPhonon);
  band.scatType.push_back(ElectronCollisionTypeOpticalPhonon);
  band.nLevels += 2;
  return true;
}

bool MediumSilicon::IntervalleyScatteringRates(
    const double rho, const double kbt, const double dtk, const double eph,
    const int kTgt, const double zTgt, const double eMinTgt,
    const int collisionType, Band& band) {
  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705
  //  - M. Lundstrom, Fundamentals of carrier transport
  //  - M. Martin et al.,
  //    Semicond. Sci. Technol. 8, 1291-1297

  // Phonon occcupation numbers
  const double nocc = 1. / (exp(eph / kbt) - 1.);
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  const double c = zTgt * c0 * dtk * dtk / eph;

  const double eMinSrc = band.eMin;
  for (int i = 0; i < band.nEnergySteps; ++i) {
    const double en = i * band.eStep;
    // Absorption
    if (en > eMinSrc && en + eph > eMinTgt) {
      band.cf[i].push_back(c * nocc * ConductionBandDOS(en + eph, kTgt));
    } else {
      band.cf[i].push_back(0.);
    }
    // Emission
    if (en > eMinSrc && en - eph > eMinTgt) {
      band.cf[i].push_back(c * (nocc + 1) * ConductionBandDOS(en - eph, kTgt));
    } else {
      band.cf[i].push_back(0.);
    }
  }
  // Absorption
  band.energyLoss.push_back(-eph);
  band.scatType.push_back(collisionType);
  // Emission
  band.energyLoss.push_back(eph);
  band.scatType.push_back(collisionType);
  band.nLevels += 2;
  return true;
}

bool MediumSilicon::IonisationRates(const std::vector<double>& p,
                                    const std::vector<double>& eth,
                                    Band& band) {
  // References:
  // - E. Cartier, M. V. Fischetti, E. A. Eklund and F. R. McFeely,
  //   Appl. Phys. Lett 62, 3339-3341
  // - DAMOCLES web page: www.research.ibm.com/DAMOCLES

  const size_t nTerms = p.size();
  if (nTerms != eth.size()) return false;
  for (int i = 0; i < band.nEnergySteps; ++i) {
    const double en = i * band.eStep;
    if (en < band.eMin) {
      band.cf[i].push_back(0.);
      continue;
    }
    double fIon = 0.;
    for (size_t j = 0; j < nTerms; ++j) {
      if (en > eth[j]) fIon += p[j] * (en - eth[j]) * (en - eth[j]);
    }
    band.cf[i].push_back(fIon);
  }

  band.energyLoss.push_back(eth[0]);
  band.scatType.push_back(ElectronCollisionTypeIonisation);
  band.nLevels += 1;
  return true;
}

bool MediumSilicon::ImpurityScatteringRates(
    const double kbt, Band& band) {

  // Density of states effective mass.
  const double md = ElectronMass * pow(band.mL * band.mT * band.mT, 1. / 3.);
  // Dielectric constant
  const double eps = GetDielectricConstant();
  // Impurity concentration
  const double impurityConcentration = m_dopingConcentration;
  if (impurityConcentration < Small) return true;

  // Screening length
  const double ls = sqrt(eps * kbt / (4 * Pi * FineStructureConstant * HbarC *
                         impurityConcentration));
  const double eb = 0.5 * HbarC * HbarC / (md * ls * ls);

  // Prefactor
  // const double c = pow(2., 2.5) * Pi * impurityConcentration *
  //                  pow(FineStructureConstant * HbarC, 2) *
  //                  SpeedOfLight / (eps * eps * sqrt(md) * eb * eb);
  // Use momentum-transfer cross-section
  const double c = impurityConcentration * Pi *
                   pow(FineStructureConstant * HbarC, 2) * SpeedOfLight /
                   (sqrt(2 * md) * eps * eps);

  for (int i = 0; i < band.nEnergySteps; ++i) {
    const double en = i * band.eStep;
    const double gamma = en * (1. + band.alpha * en);
    if (en <= band.eMin || gamma <= 0.) {
      band.cf[i].push_back(0.);
    } else {
      const double b = 4 * gamma / eb;
      band.cf[i].push_back((c / pow(gamma, 1.5)) *
                           (log(1. + b) - b / (1. + b)));
    }
  }
  band.energyLoss.push_back(0.);
  band.scatType.push_back(ElectronCollisionTypeImpurity);
  band.nLevels += 1;
  return true;
}

bool MediumSilicon::HoleScatteringRates() {
  // Reset the scattering rates
  m_vb.cfTot.assign(m_vb.nEnergySteps, 0.);
  m_vb.cf.assign(m_vb.nEnergySteps, std::vector<double>());
  m_vb.energyLoss.clear();
  m_vb.scatType.clear();
  m_vb.cfNull = 0.;

  m_vb.nLevels = 0;
  // Fill the scattering rates table
  HoleAcousticScatteringRates();
  HoleOpticalScatteringRates();
  // HoleImpurityScatteringRates();
  HoleIonisationRates();

  std::ofstream outfile;
  if (m_cfOutput) {
    outfile.open("ratesV.txt", std::ios::out);
  }

  for (int i = 0; i < m_vb.nEnergySteps; ++i) {
    // Sum up the scattering rates of all processes.
    for (int j = 0; j < m_vb.nLevels; ++j) m_vb.cfTot[i] += m_vb.cf[i][j];

    if (m_cfOutput) {
      outfile << i * m_vb.eStep << " " << m_vb.cfTot[i] << " ";
      for (int j = 0; j < m_vb.nLevels; ++j) {
        outfile << m_vb.cf[i][j] << " ";
      }
      outfile << "\n";
    }

    if (m_vb.cfTot[i] > m_vb.cfNull) {
      m_vb.cfNull = m_vb.cfTot[i];
    }

    // Make sure the total scattering rate is positive.
    if (m_vb.cfTot[i] <= 0.) {
      std::cerr << m_className << "::HoleScatteringRates:\n"
                << "    Scattering rate at " << i * m_vb.eStep << " eV <= 0.\n";
      return false;
    }
    // Normalise the rates.
    for (int j = 0; j < m_vb.nLevels; ++j) {
      m_vb.cf[i][j] /= m_vb.cfTot[i];
      if (j > 0) m_vb.cf[i][j] += m_vb.cf[i][j - 1];
    }
  }

  if (m_cfOutput) {
    outfile.close();
  }

  return true;
}

bool MediumSilicon::HoleAcousticScatteringRates() {
  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705
  //  - DAMOCLES web page: www.research.ibm.com/DAMOCLES
  //  - M. Lundstrom, Fundamentals of carrier transport

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  // Acoustic phonon intraband scattering
  // Acoustic deformation potential [eV]
  // DAMOCLES: 4.6 eV; Lundstrom: 5 eV
  constexpr double dp = 4.6;
  // Longitudinal velocity of sound [cm/ns]
  constexpr double u = 9.04e-4;
  // Prefactor for acoustic deformation potential scattering
  const double c = TwoPi * SpeedOfLight * SpeedOfLight * kbt * dp * dp /
                   (Hbar * u * u * rho);

  // Fill the scattering rate tables.
  for (int i = 0; i < m_vb.nEnergySteps; ++i) {
    const double en = Small + i * m_vb.eStep;
    m_vb.cf[i].push_back(c * ValenceBandDOS(en, 0));
  }

  // Assume that energy loss is negligible.
  m_vb.energyLoss.push_back(0.);
  m_vb.scatType.push_back(ElectronCollisionTypeAcousticPhonon);
  ++m_vb.nLevels;

  return true;
}

bool MediumSilicon::HoleOpticalScatteringRates() {
  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705
  //  - DAMOCLES web page: www.research.ibm.com/DAMOCLES
  //  - M. Lundstrom, Fundamentals of carrier transport

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  // Coupling constant [eV/cm]
  // DAMOCLES: 6.6, Lundstrom: 6.0
  constexpr double dtk = 6.6e8;
  // Phonon energy [eV]
  constexpr double eph = 63.0e-3;
  // Phonon cccupation numbers
  const double nocc = 1. / (exp(eph / kbt) - 1);
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c = c0 * dtk * dtk / eph;

  for (int i = 0; i < m_vb.nEnergySteps; ++i) {
    const double en = i * m_vb.eStep;
    // Absorption
    m_vb.cf[i].push_back(c * nocc * ValenceBandDOS(en + eph, 0));
    // Emission
    if (en > eph) {
      m_vb.cf[i].push_back(c * (nocc + 1) * ValenceBandDOS(en - eph, 0));
    } else {
      m_vb.cf[i].push_back(0.);
    }
  }

  // Absorption
  m_vb.energyLoss.push_back(-eph);
  // Emission
  m_vb.energyLoss.push_back(eph);
  m_vb.scatType.push_back(ElectronCollisionTypeOpticalPhonon);
  m_vb.scatType.push_back(ElectronCollisionTypeOpticalPhonon);

  m_vb.nLevels += 2;

  return true;
}

bool MediumSilicon::HoleIonisationRates() {
  // References:
  //  - DAMOCLES web page: www.research.ibm.com/DAMOCLES

  // Coefficients [ns-1]
  constexpr double p[2] = {2., 1.e3};
  // Threshold energies [eV]
  constexpr double eth[2] = {1.1, 1.45};
  // Exponents
  constexpr double b[2] = {6., 4.};

  double en = 0.;
  for (int i = 0; i < m_vb.nEnergySteps; ++i) {
    const double en = i * m_vb.eStep;
    double fIon = 0.;
    for (size_t j = 0; j < 2; ++j) {
      if (en > eth[j]) fIon += p[j] * pow(en - eth[j], b[j]);
    }
    m_vb.cf[i].push_back(fIon);
  }

  m_vb.energyLoss.push_back(eth[0]);
  m_vb.scatType.push_back(ElectronCollisionTypeIonisation);
  ++m_vb.nLevels;

  return true;
}

double MediumSilicon::ConductionBandDOS(const double e, const int band) {
  if (band < 0) {
    int iE = int(e / m_eStepDos);
    const int nPoints = m_fbDosConduction.size();
    if (iE >= nPoints || iE < 0) {
      return 0.;
    } else if (iE == nPoints - 1) {
      return m_fbDosConduction[nPoints - 1];
    }

    const double dos = m_fbDosConduction[iE] +
                       (m_fbDosConduction[iE + 1] - m_fbDosConduction[iE]) *
                           (e / m_eStepDos - iE);
    return dos * 1.e21;

  }
  const size_t kG = m_cb[0].nValleys + m_cb[1].nValleys;
  const auto k = m_cbIndex[band]; 
  if (k == 0) {
    // X valleys
    if (e <= 0.) return 0.;
    // Density-of-states effective mass (cube)
    const double md3 = pow(ElectronMass, 3) * m_cb[0].mL * m_cb[0].mT * m_cb[0].mT;

    if (m_fullBandDos) {
      if (e < m_cb[1].eMin) {
        return ConductionBandDOS(e, -1) / m_cb[0].nValleys;
      } else if (e < m_cb[2].eMin) {
        // Subtract the fraction of the full-band density of states
        // attributed to the L valleys.
        const double dosX =
            ConductionBandDOS(e, -1) -
            ConductionBandDOS(e, m_cb[0].nValleys) * m_cb[1].nValleys;
        return dosX / m_cb[0].nValleys;
      } else {
        // Subtract the fraction of the full-band density of states
        // attributed to the L valleys and the higher bands.
        const double dosX =
            ConductionBandDOS(e, -1) -
            ConductionBandDOS(e, m_cb[0].nValleys) * m_cb[1].nValleys -
            ConductionBandDOS(e, m_cb[0].nValleys + m_cb[1].nValleys);
        if (dosX <= 0.) return 0.;
        return dosX / m_cb[0].nValleys;
      }
    } else if (m_nonParabolic) {
      return sqrt(md3 * e * 0.5 * (1. + m_cb[0].alpha * e)) * (1. + 2 * m_cb[0].alpha * e) /
             (Pi2 * pow(HbarC, 3.));
    } else {
      return sqrt(md3 * e / 2.) / (Pi2 * pow(HbarC, 3.));
    }
  } else if (band < m_cb[0].nValleys + m_cb[1].nValleys) {
    // L valleys
    if (e <= m_cb[1].eMin) return 0.;

    // Density-of-states effective mass (cube)
    const double md3 = pow(ElectronMass, 3) * m_cb[1].mL * m_cb[1].mT * m_cb[1].mT;
    // Non-parabolicity parameter
    const double alpha = m_cb[1].alpha;

    if (m_fullBandDos) {
      // Energy up to which the non-parabolic approximation is used.
      const double ej = m_cb[1].eMin + 0.5;
      if (e <= ej) {
        return sqrt(md3 * (e - m_cb[1].eMin) * (1. + alpha * (e - m_cb[1].eMin))) *
               (1. + 2 * alpha * (e - m_cb[1].eMin)) /
               (Sqrt2 * Pi2 * pow(HbarC, 3.));
      } else {
        // Fraction of full-band density of states attributed to L valleys
        double fL = sqrt(md3 * (ej - m_cb[1].eMin) * (1. + alpha * (ej - m_cb[1].eMin))) *
                    (1. + 2 * alpha * (ej - m_cb[1].eMin)) /
                    (Sqrt2 * Pi2 * pow(HbarC, 3.));
        fL = m_cb[1].nValleys * fL / ConductionBandDOS(ej, -1);

        double dosXL = ConductionBandDOS(e, -1);
        if (e > m_cb[2].eMin) {
          dosXL -=
              ConductionBandDOS(e, m_cb[0].nValleys + m_cb[1].nValleys);
        }
        if (dosXL <= 0.) return 0.;
        return fL * dosXL / 8.;
      }
    } else if (m_nonParabolic) {
      return sqrt(md3 * (e - m_cb[1].eMin) * (1. + alpha * (e - m_cb[1].eMin))) *
             (1. + 2 * alpha * (e - m_cb[1].eMin)) / (Sqrt2 * Pi2 * pow(HbarC, 3.));
    } else {
      return sqrt(md3 * (e - m_cb[1].eMin) / 2.) / (Pi2 * pow(HbarC, 3.));
    }
  } else if (band == m_cb[0].nValleys + m_cb[1].nValleys) {
    // Higher bands
    const double ej = 2.7;
    if (m_cb[2].eMin >= ej) {
      std::cerr << m_className << "::ConductionBandDOS:\n"
                << "    Cannot determine higher band density-of-states.\n"
                << "    Program bug. Check offset energy!\n";
      return 0.;
    }
    if (e < m_cb[2].eMin) {
      return 0.;
    } else if (e < ej) {
      // Coexistence of XL and higher bands.
      const double dj = ConductionBandDOS(ej, -1);
      // Assume linear increase of density-of-states.
      return dj * (e - m_cb[2].eMin) / (ej - m_cb[2].eMin);
    } else {
      return ConductionBandDOS(e, -1);
    }
  }

  std::cerr << m_className << "::ConductionBandDOS:\n"
            << "    Band index (" << band << ") out of range.\n";
  return ElectronMass * sqrt(ElectronMass * e / 2.) / (Pi2 * pow(HbarC, 3.));
}

double MediumSilicon::ValenceBandDOS(const double e, const int band) {
  if (band <= 0) {
    // Total (full-band) density of states.
    const int nPoints = m_fbDosValence.size();
    int iE = int(e / m_eStepDos);
    if (iE >= nPoints || iE < 0) {
      return 0.;
    } else if (iE == nPoints - 1) {
      return m_fbDosValence[nPoints - 1];
    }

    const double dos =
        m_fbDosValence[iE] +
        (m_fbDosValence[iE + 1] - m_fbDosValence[iE]) * (e / m_eStepDos - iE);
    return dos * 1.e21;
  }

  std::cerr << m_className << "::ConductionBandDOS:\n"
            << "    Band index (" << band << ") out of range.\n";
  return 0.;
}

void MediumSilicon::ComputeSecondaries(const double e0, double& ee,
                                       double& eh) {
  const int nPointsValence = m_fbDosValence.size();
  const int nPointsConduction = m_fbDosConduction.size();
  const double widthValenceBand = m_eStepDos * nPointsValence;
  const double widthConductionBand = m_eStepDos * nPointsConduction;

  bool ok = false;
  while (!ok) {
    // Sample a hole energy according to the valence band DOS.
    eh = RndmUniformPos() * std::min(widthValenceBand, e0);
    int ih = std::min(int(eh / m_eStepDos), nPointsValence - 1);
    while (RndmUniform() > m_fbDosValence[ih] / m_fbDosMaxV) {
      eh = RndmUniformPos() * std::min(widthValenceBand, e0);
      ih = std::min(int(eh / m_eStepDos), nPointsValence - 1);
    }
    // Sample an electron energy according to the conduction band DOS.
    ee = RndmUniformPos() * std::min(widthConductionBand, e0);
    int ie = std::min(int(ee / m_eStepDos), nPointsConduction - 1);
    while (RndmUniform() > m_fbDosConduction[ie] / m_fbDosMaxC) {
      ee = RndmUniformPos() * std::min(widthConductionBand, e0);
      ie = std::min(int(ee / m_eStepDos), nPointsConduction - 1);
    }
    // Calculate the energy of the primary electron.
    const double ep = e0 - m_bandGap - eh - ee;
    if (ep < Small) continue;
    if (ep > 5.) return;
    // Check if the primary electron energy is consistent with the DOS.
    int ip = std::min(int(ep / m_eStepDos), nPointsConduction - 1);
    if (RndmUniform() > m_fbDosConduction[ip] / m_fbDosMaxC) continue;
    ok = true;
  }
}

void MediumSilicon::InitialiseDensityOfStates() {
  m_eStepDos = 0.1;

  m_fbDosValence = {
      {0.,      1.28083,  2.08928, 2.70763, 3.28095, 3.89162, 4.50547, 5.15043,
       5.89314, 6.72667,  7.67768, 8.82725, 10.6468, 12.7003, 13.7457, 14.0263,
       14.2731, 14.5527,  14.8808, 15.1487, 15.4486, 15.7675, 16.0519, 16.4259,
       16.7538, 17.0589,  17.3639, 17.6664, 18.0376, 18.4174, 18.2334, 16.7552,
       15.1757, 14.2853,  13.6516, 13.2525, 12.9036, 12.7203, 12.6104, 12.6881,
       13.2862, 14.0222,  14.9366, 13.5084, 9.77808, 6.15266, 3.47839, 2.60183,
       2.76747, 3.13985,  3.22524, 3.29119, 3.40868, 3.6118,  3.8464,  4.05776,
       4.3046,  4.56219,  4.81553, 5.09909, 5.37616, 5.67297, 6.04611, 6.47252,
       6.9256,  7.51254,  8.17923, 8.92351, 10.0309, 11.726,  16.2853, 18.2457,
       12.8879, 7.86019,  6.02275, 5.21777, 4.79054, 3.976,   3.11855, 2.46854,
       1.65381, 0.830278, 0.217735}};

  m_fbDosConduction = {
      {0.,      1.5114,  2.71026,  3.67114,  4.40173, 5.05025, 5.6849,  6.28358,
       6.84628, 7.43859, 8.00204,  8.80658,  9.84885, 10.9579, 12.0302, 13.2051,
       14.6948, 16.9879, 18.4492,  18.1933,  17.6747, 16.8135, 15.736,  14.4965,
       13.1193, 12.1817, 12.6109,  15.3148,  19.4936, 23.0093, 24.4106, 22.2834,
       19.521,  18.9894, 18.8015,  17.9363,  17.0252, 15.9871, 14.8486, 14.3797,
       14.2426, 14.3571, 14.7271,  14.681,   14.3827, 14.2789, 14.144,  14.1684,
       14.1418, 13.9237, 13.7558,  13.5691,  13.4567, 13.2693, 12.844,  12.4006,
       12.045,  11.7729, 11.3607,  11.14,    11.0586, 10.5475, 9.73786, 9.34423,
       9.4694,  9.58071, 9.6967,   9.84854,  10.0204, 9.82705, 9.09102, 8.30665,
       7.67306, 7.18925, 6.79675,  6.40713,  6.21687, 6.33267, 6.5223,  6.17877,
       5.48659, 4.92208, 4.44239,  4.02941,  3.5692,  3.05953, 2.6428,  2.36979,
       2.16273, 2.00627, 1.85206,  1.71265,  1.59497, 1.46681, 1.34913, 1.23951,
       1.13439, 1.03789, 0.924155, 0.834962, 0.751017}};

  auto it = std::max_element(m_fbDosValence.begin(), m_fbDosValence.end());
  m_fbDosMaxV = *it;
  it = std::max_element(m_fbDosConduction.begin(), m_fbDosConduction.end());
  m_fbDosMaxC = *it;
}
}
