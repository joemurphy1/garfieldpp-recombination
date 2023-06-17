#ifndef G_TRACK_DEGRADE_H
#define G_TRACK_DEGRADE_H

#include <array>

#include "Track.hh"

namespace Garfield {

/// Interface to Degrade.

class TrackDegrade : public Track {
 public:
  struct Electron {
    double x = 0.;
    double y = 0.;
    double z = 0.;
    double t = 0.;
    double energy = 0.;
    double dx = 0.;
    double dy = 0.;
    double dz = 0.;
  };
  struct Cluster {
    double x, y, z, t;
    double energy;
    std::vector<Electron> deltaElectrons;
    std::vector<Electron> electrons;
  };

  /// Constructor
  TrackDegrade();
  /// Destructor
  virtual ~TrackDegrade() {}

  // double GetClusterDensity() override;
  // double GetStoppingPower() override;

  bool NewTrack(const double x0, const double y0, const double z0,
                        const double t0, const double dx0, const double dy0,
                        const double dz0) override;
  bool GetCluster(double& xc, double& yc, double& zc,
                  double& tc, int& ne, double& ec, double& extra) override;
  const std::vector<Cluster>& GetClusters() const { return m_clusters; }

  void SetParticle(const std::string& particle) override;

  bool Initialise(Medium* medium, const bool verbose = false);

 protected:
  std::vector<Cluster> m_clusters;
  size_t m_cluster = 0;

  bool m_penning = true;
  bool m_bremsStrahlung = true;
  bool m_fullCascade = true;

  double m_mediumDensity = -1.;
  std::string m_mediumName = "";
  unsigned int m_nGas = 0;
 
  std::array<double, 6> m_rPenning;
  std::array<double, 6> m_dPenning;

  std::vector<Electron> TransportDeltaElectron(
      const double x0, const double y0, const double z0, const double t0,
      const double e0, const double dx, const double dy, const double dz);

  void SetupPenning(Medium* medium, std::array<double, 6> rP,
                    std::array<double, 6> dP); 
};
}

#endif
