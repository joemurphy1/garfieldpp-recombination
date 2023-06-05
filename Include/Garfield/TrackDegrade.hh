#ifndef G_TRACK_DEGRADE_H
#define G_TRACK_DEGRADE_H

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

  bool m_penning = false;
  bool m_bremsStrahlung = true;
};
}

#endif
