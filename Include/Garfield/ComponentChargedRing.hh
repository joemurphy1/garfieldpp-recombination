#ifndef G_COMPONENT_CHARGED_RING_H
#define G_COMPONENT_CHARGED_RING_H

#include <array>
#include <cmath>
#include <string>

#include "Garfield/Component.hh"

namespace Garfield {

/// Component for calculating the field of a system of charged r,z rings.

class ComponentChargedRing : public Component {
 public:
  /// Constructor
  ComponentChargedRing();
  /// Destructor
  ~ComponentChargedRing() {}

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
  
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override;

  /// Add a ring of charge N*e at x,y,z                    
  bool AddChargedRing(const double x, const double y, const double z,
                      const double N);

  /// function to set the tolerance below which two rings are considered the
  /// same.
  void SetSpacingTolerance(const double spacing_tolerance) {
    m_dSpacingTolerance = spacing_tolerance;
  }

  /// function to set the tolerance below which divergences are caught
  void SetSelfFieldTolerance(const double self_field_tolerance) {
    m_dSelfFieldTolerance = self_field_tolerance;
  }

  struct Ring {
    double z, r, charge;
    Ring(double z_, double r_, double charge_) {
      z = z_;
      r = r_;
      // This is charge / (2Pi * 4PiEpsilon0) which allows for more
      // efficient field calculation
      charge = charge_;
    }
  };

  /// Set the axis of symmetry
  void UpdateCentre(double x, double y) {
    m_centre = {x, y};
    m_bCentreSet = true;
  }

  void EnableDebugging() { m_bDebug = true; }

  void ClearActiveRings() { m_vRings.clear(); }

  std::size_t GetNumberOfRings() const { return m_vRings.size(); }

  const std::vector<Ring>& GetRings() const { return m_vRings; }

 private:
  /// Active area.
  std::array<double, 3> m_xmin = {{0., 0., 0.}};
  std::array<double, 3> m_xmax = {{0., 0., 0.}};
  /// Did we specify the active area explicitly?
  bool m_hasArea = false;
  /// Medium in the active area.
  Medium* m_medium = nullptr;

  void Reset() override;
  void UpdatePeriodicity() override;
  
  /// Gets elliptic integrals via list
  void GetEllipticIntegrals(double x, double& K, double& E) const;

  bool InArea(const double x, const double y, const double z) {
    if (x < m_xmin[0] || x > m_xmax[0] || y < m_xmin[1] || y > m_xmax[1] ||
        z < m_xmin[2] || z > m_xmax[2]) {
      return false;
    }
    return true;
  }

  enum class Elliptic : std::size_t { X, K, E };
  static const constexpr std::size_t elliptic_size{29981};
  static const std::array<std::array<double, 3>, elliptic_size> m_elliptic;

  bool m_bDebug = false;

  /// centre of cylindrical symmetry (x,y)
  std::array<double, 2> m_centre = {0., 0.};

  std::vector<Ring> m_vRings;

  double m_dSelfFieldTolerance = 0.00001;  // this is a decent value with
                                           // minimal 'strangeness' in the field
  double m_dSpacingTolerance = 2.1 * m_dSelfFieldTolerance;

  bool m_bCentreSet = false;

  static void GetCoulombBallField(const Ring& ring, 
                                  const double r, const double z,
                                  double& eFieldZ, double& eFieldR);
  void GetChargedRingField(const Ring& ring, double r, double z,
                           double& eFieldZ, double& eFieldR) const;

};
}  // namespace Garfield
#endif
