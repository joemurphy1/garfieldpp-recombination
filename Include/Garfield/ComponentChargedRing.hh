#ifndef G_COMPONENT_CONSTANT_H
#define G_COMPONENT_CONSTANT_H

#include <array>
#include <string>
#include <cmath>

#include "Garfield/Component.hh"

namespace Garfield {

/// Component for calculating the field of a system of charged r,z rings.

class ComponentChargedRing : public Component {
 public:
  /// Constructor
  ComponentChargedRing();
  /// Destructor
  ~ComponentChargedRing() {}

  /// Set the components of the electric field [V / cm].
  void SetElectricField(const double ex, const double ey, const double ez);
  /// Specify the potential at a given point.
  void SetPotential(const double x, const double y, const double z,
                    const double v = 0.);

  /// Set the components of the weighting field [1 / cm].
  void SetWeightingField(const double wx, const double wy, const double wz,
                         const std::string label);
  /// Specify the weighting potential at a given point.
  void SetWeightingPotential(const double x, const double y, const double z,
                             const double v = 0.);

  /// Set the limits of the active area explicitly
  /// (instead of using a Geometry object).
  void SetArea(const double xmin, const double ymin, const double zmin,
               const double xmax, const double ymax, const double zmax);
  /// Remove the explicit limits of the active area.
  void UnsetArea();
  /// Set the medium in the active area.
  void SetMedium(Medium* medium) { m_medium = medium; }

  Medium* GetMedium(const double x, const double y, const double z) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;

  using Component::ElectricField;
  bool GetVoltageRange(double& vmin, double& vmax) override;
  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override;

  bool AddChargedRing(const double x, const double y, const double z, const int N);

  // function to set the tolerance below which two rings are considered the same.
  void SetSpacingTolerance(const double spacing_tolerance){
    m_dSpacingTolerance = spacing_tolerance;
  }

  // function to set the tolerance below which the self field is rejected
  void SetSelfFieldTolerance(const double self_field_tolerance){
    m_dSelfFieldTolerance = self_field_tolerance;
  }

  struct Ring{
    double z,r,charge; // charge in coulombs
    Ring(double z_,double r_,double charge_){
        z = z_;
        r = r_;
        charge = charge_;
    }
  };

  void UpdateCentre(double x, double y){
    m_centre = {x,y};
    m_bCentreSet = true;
  }

  void EnableDebugging(){
    m_bDebug = true;
  }

  void ClearActiveRings(){
    m_vRings.clear();
  }

  void GetNumberOfRings(int & n_rings){
    n_rings = m_vRings.size();
  }

 private:

  // Active area.
  std::array<double, 3> m_xmin = {{0., 0., 0.}};
  std::array<double, 3> m_xmax = {{0., 0., 0.}};
  // Did we specify the active area explicitly?
  bool m_hasArea = false;
  // Medium in the active area.
  Medium* m_medium = nullptr;

  void Reset() override;
  void UpdatePeriodicity() override;

  bool InArea(const double x, const double y, const double z) {
    if (x < m_xmin[0] || x > m_xmax[0] || y < m_xmin[1] || y > m_xmax[1] ||
        z < m_xmin[2] || z > m_xmax[2]) {
      return false;
    }
    return true;
  }


  // Import elliptic integral values
  void ImportEllipticIntegralValues(const std::string &filename);

  bool m_bImportElliptic = false;

  enum class Elliptic : std::size_t
  {
    X,
    K,
    E
  };
  static const constexpr std::size_t elliptic_size{29981};
  static const std::array<std::array<double,3>, elliptic_size> m_elliptic;

  // Gets elliptic integrals via list
  void GetEllipticIntegrals(double x, double &K, double &E);

  bool m_bDebug = false;

  // centre of cylindrical symmetry (x,y)
  std::array<double,3> m_centre = {0.,0.};

  std::vector<Ring> m_vRings;

  double m_dSelfFieldTolerance = 0.00001; //this is a decent value with minimal 'strangeness' in the field
  double m_dSpacingTolerance = 2.1*m_dSelfFieldTolerance; // > twice SFT to avoid divergence catching converting balls to rings
  
  bool m_bCentreSet = false;

  void GetCoulombBallField(const Ring & ring, const double r, const double z, double & eFieldZ, double & eFieldR);

  void GetChargedRingField(const Ring & ring, const double r, const double z, double & eFieldZ, double & eFieldR);

  void GetCartesianLocalField(double &eFieldR, double &eFieldX, double &eFieldY, double x, double y){
    // get cartesian field
    double cosphi = std::cos(std::atan2(y,x));
    double sinphi = std::sin(std::atan2(y,x));
    eFieldX = eFieldR*cosphi;
    eFieldY = eFieldR*sinphi;
  }

};
}  // namespace Garfield
#endif
