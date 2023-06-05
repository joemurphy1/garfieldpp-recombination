#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackDegrade.hh"
#include "Garfield/DegradeInterface.hh"

namespace {

Garfield::TrackDegrade::Electron MakeElectron(
    const double energy, const double x, const double y, const double z, 
    const double t, const double dx, const double dy, const double dz) {

  Garfield::TrackDegrade::Electron electron;
  electron.energy = energy;
  electron.x = x;
  electron.y = y;
  electron.z = z;
  electron.t = t;
  electron.dx = dx;
  electron.dy = dy;
  electron.dz = dz;
  return electron;
}

}

namespace Garfield {

TrackDegrade::TrackDegrade() : Track("Degrade") {

  m_q = -1;
  m_spin = 1;
  m_mass = ElectronMass;
  m_isElectron = true;
  SetBetaGamma(3.);
  m_particleName = "electron";
}

bool TrackDegrade::NewTrack(const double x0, const double y0, const double z0,
                            const double t0, const double dx0, const double dy0,
                            const double dz0) {

  m_clusters.clear();
  m_cluster = 0;
  // Make sure the sensor is defined.
  if (!m_sensor) {
    std::cerr << m_className << "::NewTrack: Sensor is not defined.\n";
    return false;
  }

  // Make sure we are inside a medium.
  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::NewTrack: No medium at initial position.\n";
    return false;
  }
  if (!Initialise(medium, m_debug)) return false;

  double xp = x0;
  double yp = y0;
  double zp = z0;
  double tp = t0;

  // Normalise the direction.
  double dxp = dx0;
  double dyp = dy0;
  double dzp = dz0;
  const double dp = sqrt(dxp * dxp + dyp * dyp + dzp * dzp);
  if (dp < Small) {
    // Choose a random direction.
    RndmDirection(dxp, dyp, dzp);
  } else {
    const double scale = 1. / dp;
    dxp *= scale;
    dyp *= scale;
    dzp *= scale;
  }
  // Get the total collision rate.
  int64_t ie = 20000;
  double tcf = 0.;
  Degrade::gettcf(&ie, &tcf);
  // Convert to ns.
  tcf *= 1.e3;
  double beta = GetBeta();
  double gamma = GetGamma();
  double ep = GetKineticEnergy();
  if (m_debug) {
    std::cout << m_className << "::NewTrack:\n"
              << "    Particle energy: " << ep << " eV\n"
              << "    Collision rate: " << tcf << " / ns\n"
              << "    Collisions / cm: " 
              << 1. / (SpeedOfLight * beta / tcf) << "\n";
  }
  bool ok = true;
  while (ok) {
    // Draw a time step.
    double dt = -log(RndmUniformPos()) / tcf;
    const double step = SpeedOfLight * beta * dt;
    xp += dxp * step;
    yp += dyp * step;
    zp += dzp * step;
    tp += dt;
    if (!m_sensor->IsInside(xp, yp, zp)) break;
    // Reset the scattering angle.
    double cthetap = -99.;
    double sthetap = -99.;
    // Determine the type of interaction.
    double r1 = RndmUniform();
    int64_t izbr = 0;
    double rgas = 0.;
    double ein = 0.;
    int64_t ia = 0;
    double wpl = 0.;
    int64_t index = 0;
    double an = 0.;
    double ps = 0.;
    double wklm = 0.;
    int64_t nc0 = 0;
    double ec0 = 0.;
    int64_t ng2 = 0;
    double eg2 = 0.;
    int64_t ng1 = 0;
    double eg1 = 0.;
    double dstfl = 0.;
    Degrade::getlevel(&ie, &r1, &izbr, &rgas, &ein, &ia, &wpl, &index, 
                      &an, &ps, &wklm, &nc0, &ec0, &ng1, &eg1, &ng2, &eg2, 
                      &dstfl); 
    // Bremsstrahlung?
    if (izbr != 0 && m_bremsStrahlung) {
      double eout = 0., egamma = 0.;
      double dxe = 0., dye = 0., dze = 0.; 
      double dxg = 0., dyg = 0., dzg = 0.;
      Degrade::brems(&izbr, &ep, &dxp, &dyp, &dzp, &eout, &dxe, &dye, &dze,
                     &egamma, &dxg, &dyg, &dzg);
      ep = eout;
      dxp = dxe;
      dyp = dye;
      dzp = dze;
      continue;
    }

    if (ia ==  2 || ia ==  7 || ia == 12 || ia == 17 || 
        ia == 22 || ia == 27) {
      Cluster cluster;
      cluster.x = xp;
      cluster.y = yp;
      cluster.z = zp;
      cluster.t = tp;
      // Ionisation
      double esec = wpl * tan(RndmUniform() * atan((ep - ein) / (2. * wpl)));
      esec = wpl * pow(esec / wpl, 0.9524);
      // TODO: while esec > ECUT?
      // Calculate primary scattering angle.
      if (index == 1) {
        cthetap = 1. - RndmUniform() * an;
        if (RndmUniform() > ps) cthetap = -cthetap;
      } else if (index == 2) {
        const double r3 = RndmUniform();
        cthetap = 1. - (2. * r3 * (1. - ps) / (1. + ps * (1. - 2. * r3)));
      } else { 
        cthetap = 1. - 2. * RndmUniform();
      }
      sthetap = sin(acos(cthetap));
      const double gammas = (ElectronMass + esec) / ElectronMass;
      // Calculate secondary recoil angle from free kinematics.
      const double sthetas = std::min(sthetap * sqrt(ep / esec) * 
                                      gamma / gammas, 1.);
      double thetas = asin(sthetas);
      double phis = TwoPi * RndmUniform();
      // Calculate new direction cosines from initial values and 
      // scattering angles.
      double dxs = 0., dys = 0., dzs = 0.;
      Degrade::drcos(&dxp, &dyp, &dzp, &thetas, &phis, &dxs, &dys, &dzs);
      // Add the secondary to the list.
      cluster.electrons.emplace_back(
        MakeElectron(esec, xp, yp, zp, tp, dxs, dys, dzs)); 
      // Calculate possible shell emissions (Auger or fluorescence).
      if (wklm > 0.0 && RndmUniform() < wklm) {
        // Auger emission and fluorescence.
        if (ng2 > 0) {
          const double eavg = eg2 / ng2;
          for (int64_t j = 0; j < ng2; ++j) {
            // Random emission angle.
            const double ctheta = 1. - 2. * RndmUniform();
            const double stheta = sin(acos(ctheta));
            const double phi = TwoPi * RndmUniform();
            const double dx = cos(phi) * stheta;
            const double dy = sin(phi) * stheta;
            const double dz = ctheta;
            cluster.electrons.emplace_back(
              MakeElectron(eavg, xp, yp, zp, tp, dx, dy, dz)); 
          }
        }
        if (ng1 > 0) {
          const double eavg = eg1 / ng1;
          // Fluorescence absorption distance.
          double dfl = -log(RndmUniformPos()) * dstfl;
          for (int64_t j = 0; j < ng1; ++j) {
            // Random emission angle.
            const double ctheta = 1. - 2. * RndmUniform();
            const double stheta = sin(acos(ctheta));
            const double phi = TwoPi * RndmUniform();
            const double dx = cos(phi) * stheta;
            const double dy = sin(phi) * stheta;
            const double dz = ctheta;
            const double cthetafl = 1. - 2. * RndmUniform();
            const double sthetafl = sin(acos(cthetafl));
            const double phifl = TwoPi * RndmUniform();
            const double xs = dfl * sthetafl * cos(phifl);
            const double ys = dfl * sthetafl * sin(phifl);
            const double zs = dfl * cthetafl; 
            const double ts = dfl / SpeedOfLight;
            cluster.electrons.emplace_back(
              MakeElectron(eavg, xs, ys, zs, ts, dx, dy, dz)); 
          }
        }
      } else {
        // Auger emission without fluorescence.
        const double eavg = ec0 / nc0;
        for (int64_t j = 0; j < nc0; ++j) {
          // Random emission angle.
          const double ctheta = 1. - 2. * RndmUniform();
          const double stheta = sin(acos(ctheta));
          const double phi = TwoPi * RndmUniform();
          const double dx = cos(phi) * stheta;
          const double dy = sin(phi) * stheta;
          const double dz = ctheta;
          cluster.electrons.emplace_back(
            MakeElectron(eavg, xp, yp, zp, tp, dx, dy, dz)); 
        }
      }
      m_clusters.push_back(std::move(cluster));
    } else if (ia ==  4 || ia ==  9 || ia == 14 || ia == 19 || 
               ia == 24 || ia == 29) {
      // Excitation
      // TODO: PENFRA(1,I), PENFRA(2,I), PENFRA(3,I)
      constexpr double penfra1 = 0.;
      constexpr double penfra2 = 0.;
      constexpr double penfra3 = 0.;
      if (m_penning && penfra1 > 0. && RndmUniform() < penfra1) {
        // Penning transfer
        Cluster cluster;
        cluster.x = xp;
        cluster.y = yp;
        cluster.z = zp;
        cluster.t = tp;
        // Random emission angle.
        const double ctheta = 1. - 2. * RndmUniform();
        const double stheta = sin(acos(ctheta));
        const double phi = TwoPi * RndmUniform();
        const double dx = cos(phi) * stheta;
        const double dy = sin(phi) * stheta;
        const double dz = ctheta;
        // Penning transfer distance.
        double asign = 1.;
        if (RndmUniform() < 0.5) asign = -asign;
        // TODO: why penfra2 * 1.D-6?
        const double xs = xp - log(RndmUniformPos()) * penfra2 * asign;
        if (RndmUniform() < 0.5) asign = -asign;
        const double ys = yp - log(RndmUniformPos()) * penfra2 * asign;
        if (RndmUniform() < 0.5) asign = -asign;
        const double zs = zp - log(RndmUniformPos()) * penfra2 * asign;
        const double ts = tp - log(RndmUniformPos()) * penfra3;
        // Fix Penning electron energy to 4 eV.
        cluster.electrons.emplace_back(
          MakeElectron(4., xs, ys, zs, ts, dx, dy, dz));
        m_clusters.push_back(std::move(cluster)); 
      }
    }
    double s1 = 1. + gamma * (rgas - 1.);
    double s2 = (s1 * s1) / (s1 - 1.); 
    if (cthetap < -1.) {
      if (index == 1) {
        // Anisotropic scattering
        cthetap = 1. - RndmUniform() * an;          
        if (RndmUniform() > ps) cthetap = -cthetap;
      } else if (index == 2) {
        // Anisotropic scattering
        const double r3 = RndmUniform();
        cthetap = 1. - (2. * r3 * (1. - ps) / (1. + ps * (1. - 2. * r3))); 
      } else { 
        // Isotropic scattering
        cthetap = 1. - 2. * RndmUniform();
      }
      sthetap = sin(acos(cthetap));
    }
    double phi0 = TwoPi * RndmUniform();
    double sphi0 = sin(phi0);
    double cphi0 = cos(phi0);
    if (ep < ein) ein = 0.;
    double arg1 = std::max(1. - s1 * ein / ep, 1.e-20);
    double d = 1. - cthetap * sqrt(arg1);
    double e1 = std::max(ep * (1. - ein / (s1 * ep) - 2. * d / s2), 1.e-20);
    double q = std::min(sqrt((ep / e1) * arg1) / s1, 1.); 
    double theta = asin(q * sthetap);
    double ctheta = cos(theta);
    if (cthetap < 0.) {
      double u = (s1 - 1.) * (s1 - 1.) / arg1;
      if (cthetap * cthetap > u) ctheta = -ctheta;
    }
    double stheta = sin(theta);
    double dx1 = dxp;
    double dy1 = dyp;
    double dz1 = std::min(dzp, 1.);
    double argz = sqrt(dx1 * dx1 + dy1 * dy1);
    if (argz == 0.) {
      dzp = ctheta;
      dxp = cphi0 * stheta;
      dyp = sphi0 * stheta;
    } else {
      dzp = dz1 * ctheta + argz * stheta * sphi0;
      dyp = dy1 * ctheta + (stheta/argz) * (dx1 * cphi0 - dy1 * dz1 * sphi0);
      dxp = dx1 * ctheta - (stheta/argz) * (dy1 * cphi0 + dx1 * dz1 * sphi0);
    }
    ep = e1;
  } 
  return true;
}

bool TrackDegrade::GetCluster(double& xc, double& yc, double& zc, double& tc,
                              int& ne, double& ec, double& extra) {
  xc = yc = zc = tc = ec = extra = 0.;
  ne = 0;
  if (m_clusters.empty() || m_cluster >= m_clusters.size()) return false;
  const auto& cluster = m_clusters[m_cluster];
  xc = cluster.x;
  yc = cluster.y;
  zc = cluster.z;
  tc = cluster.t;
  ec = cluster.energy;
  ne = 1; 
  ++m_cluster;
  return true;
}

bool TrackDegrade::Initialise(Medium* medium, const bool verbose) {

  if (!medium) {
    std::cerr << m_className << "::Initialise: Null pointer.\n";
    return false;
  }
  if (!medium->IsGas()) {
    std::cerr << m_className << "::Initialise: Medium " 
              << medium->GetName() << " is not a gas.\n";
    return false;
  }

  // Get temperature and pressure.
  double pressure = medium->GetPressure();
  double temperature = medium->GetTemperature() - ZeroCelsius;
  // Get the gas composition.
  std::array<int64_t, 6> ngas;
  std::array<double, 6> frac;
  const unsigned int nComponents = medium->GetNumberOfComponents();
  if (nComponents < 1) {
    std::cerr << m_className << "::Initialise: Invalid gas mixture.\n";
    return false;
  }
  for (unsigned int i = 0; i < nComponents; ++i) {
    std::string name;
    double f;
    medium->GetComponent(i, name, f);
    ngas[i] = MediumMagboltz::GetGasNumberMagboltz(name);
    frac[i] = 100. * f;
  }

  int64_t ng = nComponents;
  int64_t ne = 1;
  int64_t mip = 1;
  int64_t idvec = 1;
  int64_t iseed = 0;
  double e0 = GetKineticEnergy();
  double et = 2.;
  double ec = 10.;
  double etot = 100.;
  double btot = 0.;
  double bang = 90.;
  int64_t jcmp = 1;
  int64_t jray = 1;
  int64_t jpap = 1;
  int64_t jbrm = m_bremsStrahlung ? 1 : 0;
  int64_t jecasc = 1;
  int64_t iverb = verbose ? 1 : 0;
  Degrade::deginit(&ng, &ne, &mip, &idvec, &iseed, &e0, &et, &ec,
                   &ngas[0], &ngas[1], &ngas[2], &ngas[3], &ngas[4], &ngas[5],
                   &frac[0], &frac[1], &frac[2], &frac[3], &frac[4], &frac[5],
                   &temperature, &pressure, &etot, &btot, &bang,
                   &jcmp, &jray, &jpap, &jbrm, &jecasc, &iverb);
  return true;
}

void TrackDegrade::SetParticle(const std::string& particle) {
  if (particle != "electron" && particle != "e" && particle != "e-") {
    std::cerr << m_className << "::SetParticle: Only electrons are allowed.\n";
  }
}

}
