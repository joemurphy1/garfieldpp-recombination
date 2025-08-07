#include "Garfield/ComponentChargedRing.hh"

#include <array>
#include <iostream>
#include <numeric>
#include <fstream>
#include <string>
#include <sstream>

#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Medium.hh"

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Utilities.hh"

namespace Garfield {

ComponentChargedRing::ComponentChargedRing() : Component("Charged ring") {}

Medium* ComponentChargedRing::GetMedium(const double x, const double y,
                                     const double z) {
  if (!m_hasArea) return Component::GetMedium(x, y, z);
  return InArea(x, y, z) ? m_medium : nullptr;
}

void ComponentChargedRing::ElectricField(const double x, const double y,
                                      const double z, double& ex, double& ey,
                                      double& ez, Medium*& m, int& status) {
  m = GetMedium(x, y, z);
  if (!m) {
    // No medium at this point.
    status = -6;
    return;
  }

  if (m->IsDriftable()) {
    status = 0;
  } else {
    status = -5;
  }
 
  if (!m_bCentreSet) {
    std::cerr << m_className
              << "::ElectricField: Centre not set.\n";
    return;
  }

  // assume cylindrical axis is y rather than z
  double r = std::sqrt((x-m_centre[0])*(x-m_centre[0]) + (z-m_centre[1])*(z-m_centre[1]));
  double eFieldR = 0;
  double eY_temp;
  double eR_temp;
  ey = 0;
  

  for (ComponentChargedRing::Ring & ring:m_vRings){
    GetChargedRingField(ring, r, y, eY_temp,eR_temp);
    eFieldR += eR_temp;
    ey += eY_temp;
  }
  GetCartesianLocalField(eFieldR,ex,ez,(x-m_centre[0]),(z-m_centre[1]));
  
  
}

void ComponentChargedRing::ElectricField(const double x, const double y,
                                      const double z, double& ex, double& ey,
                                      double& ez, double& v, Medium*& m,
                                      int& status) {
  m = GetMedium(x, y, z);
  if (!m) {
    if (m_bDebug) {
      std::cout << m_className << "::ElectricField: No medium at (" << x << ", "
                << y << ", " << z << ").\n";
    }
    status = -6;
    return;
  }

  if (m->IsDriftable()) {
    status = 0;
  } else {
    status = -5;
  }

  ElectricField(x,y,z,ex,ey,ez,m,status);
}

bool ComponentChargedRing::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                       double& xmax, double& ymax,
                                       double& zmax) {
  if (!m_hasArea) {
    return Component::GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
  }
  xmin = m_xmin[0];
  ymin = m_xmin[1];
  zmin = m_xmin[2];
  xmax = m_xmax[0];
  ymax = m_xmax[1];
  zmax = m_xmax[2];
  return true;
}


void ComponentChargedRing::SetArea(const double xmin, const double ymin,
                                const double zmin, const double xmax,
                                const double ymax, const double zmax) {
  m_xmin[0] = std::min(xmin, xmax);
  m_xmin[1] = std::min(ymin, ymax);
  m_xmin[2] = std::min(zmin, zmax);
  m_xmax[0] = std::max(xmin, xmax);
  m_xmax[1] = std::max(ymin, ymax);
  m_xmax[2] = std::max(zmin, zmax);
  m_hasArea = true;
}

void ComponentChargedRing::UnsetArea() {
  m_xmin.fill(0.);
  m_xmax.fill(0.);
  m_hasArea = false;
}

void ComponentChargedRing::Reset() {
  UnsetArea();
  m_medium = nullptr;
  m_bCentreSet = false;
  ClearActiveRings();
  m_bToleranceSet = false;
}



void ComponentChargedRing::ImportEllipticIntegralValues(const std::string &filename) {
  // reads values of the elliptic functions
  // it has a very special form to find the values rapidly (not given for a
  // completely non uniform grid) from x = 0 to 10 it is in steps of 1e-3. From
  // 10 to 1e4 in steps of 1. Then in steps of 1000 until 1e7.

  m_vEElliptic.reserve(20000);
  m_vKElliptic.reserve(20000);
  m_vXElliptic.reserve(20000);

  m_vXElliptic.resize(0);
  m_vKElliptic.resize(0);
  m_vEElliptic.resize(0);

  std::ifstream ellipticStream(filename);

  if (!ellipticStream) {
    std::cerr << m_className
              << "::ImportEllipticIntegralValues: Could not open file.\n";
  }

  for (std::string line; std::getline(ellipticStream, line);) {
    std::istringstream iss(line);
    double value = 0.;
    iss >> value;
    m_vXElliptic.push_back(value);
    iss >> value;
    m_vKElliptic.push_back(value);
    iss >> value;
    m_vEElliptic.push_back(value);
  }

  ellipticStream.close();
  m_bImportElliptic = true;
}

void ComponentChargedRing::GetEllipticIntegrals(double x, double &K,
                                                    double &E) {
  // from x = 0 to 10 it is in steps of 1e-3. From 10 to 1e4 in steps of 1. Then
  // in steps of 1000 until 1e7.

  

  if (!m_bImportElliptic) {
    // Import the elliptic integral values
    const std::string path = std::getenv("GARFIELD_INSTALL");
    std::string filename = path + "/share/Garfield/Data/elliptic_integrals.txt";
    ImportEllipticIntegralValues(filename);
  }

  int arg;
  double invStep;
  if (-x < 1.e1) {
    invStep = 1000.;
    arg = (int)(-x * invStep);
  } else if (-x < 1.e4) {
    invStep = 1.;
    arg = (int)(-x - 10) + 10000;
  } else if (-x < 1.e7) {
    invStep = 0.001;
    arg = (int)((-x - 1.e4) * invStep) + 19990;
  } else {
    // not included in list.
    if (m_bDebug)
      std::cerr << m_className
                << "::GetEllipticIntegrals: Value not included in list.\n";
    K = m_vKElliptic.back();
    E = m_vEElliptic.back();
    return;
  }

  // Linear interpolation:
  const double f = (-x - m_vXElliptic.at(arg)) * invStep;
  K = (1. - f) * m_vKElliptic.at(arg) + f * m_vKElliptic.at(arg + 1);
  E = (1. - f) * m_vEElliptic.at(arg) + f * m_vEElliptic.at(arg + 1);
}

bool ComponentChargedRing::AddChargedRing(const double x, const double y, const double z, const int N){
    if (!m_bToleranceSet){
        std::cerr << m_className << "::AddChargedRing: Tolerance not set.\n";
        return false;
    }

    // assume cylindrical axis is y rather than z
    double r = std::sqrt((x-m_centre[0])*(x-m_centre[0]) + (z-m_centre[1])*(z-m_centre[1]));

    ComponentChargedRing::Ring ring(y,r,N*ElementaryCharge);
    bool in_list = false;
    for (ComponentChargedRing::Ring & existing_ring:m_vRings){
        if (std::abs(existing_ring.z - ring.z) < m_dSpacingTolerance && std::abs(existing_ring.r - ring.r) < m_dSpacingTolerance){
            existing_ring.charge += ring.charge;
            in_list = true;      
        }
    }
    if (!in_list) {
      m_vRings.push_back(ring);
      if (m_bDebug) std::cout << m_className << "::AddChargedRing: Added ring of charge " << N << " at z = " << ring.z << ".\n";
    }
    return true;  
}

void ComponentChargedRing::GetChargedRingField(const ComponentChargedRing::Ring & ring, const double r, const double z, double & eFieldZ, double & eFieldR){
    

    // reject if the field is being calculated on the ring
    double self_field_tolerance = 0.00001;
    if (std::abs(r-ring.r) < self_field_tolerance && std::abs(z-ring.z) < self_field_tolerance) {
        eFieldZ = 0;
        eFieldR = 0;
        return;
    }

    if (ring.r < m_dSpacingTolerance) {
        GetCoulombBallField(ring, r, z, eFieldZ, eFieldR);
        return;
    }

    double dz = z - ring.z;  //< I double-checked that's the right sign

    // parameters (see Lippmann Diss.)
    const double a2 = (r + ring.r) * (r + ring.r) + dz * dz;
    const double b2 = (r - ring.r) * (r - ring.r) + dz * dz;
    const double b = std::sqrt(b2);
    const double c2 = r * r - ring.r * ring.r - dz * dz;
    // parameter for elliptic integrals
    const double x =
        -4 * r * ring.r / b2;  //< x < 0, i.e. never near x = 1 (singularity)

    // calculation of elliptic integrals and fields (up to prefactor)
    double EllE, EllK;
    GetEllipticIntegrals(x, EllK, EllE);
    eFieldZ = EllE * 4. * dz / (a2 * b);
    eFieldR = c2 * EllE + a2 * EllK;
    // if ri = 0?
    if (r < Small) {
        eFieldR = 0;
    } else {
        eFieldR *= 2. / (r * a2 * b);
    }
 
    eFieldR *= ring.charge/(TwoPi * FourPiEpsilon0);
    eFieldZ *= ring.charge/(TwoPi * FourPiEpsilon0);

}

void ComponentChargedRing::GetCoulombBallField(const ComponentChargedRing::Ring & ring, const double r, const double z, double & eFieldZ, double & eFieldR){
    const double d = std::sqrt((z - ring.z) * (z - ring.z) + r * r);
    const double f = TwoPi / (d * d * d);
    eFieldR = f * r;
    eFieldZ = f * (z - ring.z);   

    eFieldR *= ring.charge/(TwoPi * FourPiEpsilon0);
    eFieldZ *= ring.charge/(TwoPi * FourPiEpsilon0);
}

bool ComponentChargedRing::GetVoltageRange(double& vmin, double& vmax) {
  if (m_bDebug) std::cout << "GetVoltageRange not implemented.\n";
}

void ComponentChargedRing::UpdatePeriodicity() {
  if (m_bDebug) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Periodicities are not supported.\n";
  }
}

double ComponentChargedRing::WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) {
  if (m_bDebug) std::cout << "WeightingPotential not implemented.\n";
}


void ComponentChargedRing::WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label){
  if (m_bDebug) std::cout << "WeightingField not implemented.\n";
}

}  // namespace Garfield
