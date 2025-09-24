#ifndef GARFIELD_AVALANCHEGRIDSPACECHARGE_HH
#define GARFIELD_AVALANCHEGRIDSPACECHARGE_HH

#include <array>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Garfield/ComponentChargedRing.hh"

namespace Garfield {
class Sensor;
class AvalancheMicroscopic;
class ComponentParallelPlate;

/// Propagates avalanches with the 2d (axi-symmetric) space-charge routine from
/// Lippmann, Riegler (2004) in uniform background fields. Different options to
/// calculate space-charge-fields can be chosen.
class AvalancheGridSpaceCharge {
 public:
  /// Default constructor
  AvalancheGridSpaceCharge();
  /// Constructor
  explicit AvalancheGridSpaceCharge(Sensor *sensor);
  explicit AvalancheGridSpaceCharge(std::nullptr_t) = delete;
  /// Destructor
  ~AvalancheGridSpaceCharge() = default;

  /// Set the sensor (and determine if it includes a parallel-plate component).
  void SetSensor(Sensor *sensor);

  /// Reset the charges.
  void Reset();

  /// Enable/disable debugging (log messages) (default off)
  void EnableDebugging(const bool option = true) { m_bDebug = option; }

  /// Enable sticky anode (default on)
  void EnableStickyAnode(const bool option = true) { m_bStick = option; }

  /// Enable to use TOF swarm parameters (default on)
  void EnableTOF(const bool option = true) { m_bUseTOF = option; }

  /// Enable diffusion (default off)
  void EnableDiffusion(const bool option = true) { m_bDiffusion = option; }

  /// Enable space charge calculations (default on)
  void EnableSpaceChargeEffect(const bool option = true) {
    m_bSpaceCharge = option;
    if (option) {
      if (!m_isgridset) {
        throw std::runtime_error(m_className +
                                 "::EnableSpaceChargeEffect: use Set2dGrid() "
                                 "before enabling space charge.");
      }
      SetRingSystems();
    }
  }

  /// Enable adaptive time stepping (default on)
  void EnableAdaptiveTimeStepping(const bool option = true) {
    m_bAdaptiveTime = option;
  }

  /// Enable Monte Carlo (gain/diffusion) up to a certain total number of
  /// electrons e.g. (1e5) (default on)
  ///   if disabled: Mean values are considered instead of the statistical
  ///   processes
  void EnableMC(const bool option = true) { m_bMC = option; }

  /// Set electron saturation value applied for each gas gap if space charge
  /// effect is turned off (default 1e8)
  void SetNCrit(const long NCrit = 1e8) { m_lNCrit = NCrit; }

  /// Sets the method for calculating the space charge field.
  ///   Coulomb: free field approximation.
  ///   Mirror: symmetric three-layer single-gap RPC
  ///           (metal - resistive layer - gas - resistive layer - metal).
  void SetFieldCalculation(const std::string &option = "coulomb",
                           const int nof_approx = 1);

  /// Set the streamer-inception criterion constant K in the interval (0,
  /// &infin;) s.t. 1 = 100%
  void SetK(float option = 0.95) { m_fStreamerK = option; }

  /// Stop the avalanche if K % field is reached
  void SetStopAtK(bool option = true) { m_bStopAtK = option; }

  /**
   *
   * @param zmin coordinate along direction of background field
   * @param zmax should more or less be at the end of the gap
   * @param zsteps number of steps
   * @param rmax coordinate perpendicular to background field
   * @param rsteps number of steps
   */
  void Set2dGrid(double zmin, double zmax, int zsteps, double rmax, int rsteps);

  /// Import electron (no ions) data from AvalancheMicroscopic class to
  /// axi-symmetric grid
  void AddElectrons(AvalancheMicroscopic *avmc);

  /// Set n electrons onto the grid
  void AddElectron(double x, double y, double z, double t = 0, int n = 1);

  /// After calling AddElectron, add more electrons on the same
  /// transversal line (y-freedom).
  void AddExtraElectron(double y, int n = 1);

  /// Starts the simulation with the imported electrons for a time step dt
  /// (dt = -1 until there are no electrons left in the gap).
  void StartGridAvalanche(double dtime = -1);

  /// Returns the total positive charge in the gap's
  long GetAvalancheSize() const { return m_nTotPosIons; }

  /// Return current mean distance of the electrons on the grid
  double GetMeanDistance();

  /// Returns if 100 * K % background charge has been reached
  [[nodiscard]] long ReachedKPercent() const {
    if (m_bFieldK)
      return m_lElectronsK;
    else
      return -1;
  }

  /// Returns the total electron number evolution
  [[nodiscard]] const std::vector<std::pair<double, long>> &
  GetElectronEvolution() const {
    return m_vNElectronEvolution;
  }

  /// Export the current grid to a txt file (electron, ions numbers and field
  /// magnitude)
  ///  Take care: only where electrons are located is the space-charge field
  ///  evaluated!
  void ExportGrid(const std::string &filename);

 private:
  struct GridNode {
    long nE{0};     ///< number of electrons
    double nP{0.};  ///< number of positive ions (smeared values allowed)
    double nN{0.};  ///< number of negative ions (smeared values allowed)
    // Electrons and ions at the next step in time:
    long nEHolder{0};     ///< at t+dt
    double nPHolder{0.};  ///< at t+dt
    double nNHolder{0.};  ///< at t+dt

    double townsend{0.};    ///< Townsend coefficient [1/cm]
    double attachment{0.};  ///< Attachment coefficient [1/cm]
    /// Magnitude of the drift velocity [cm/ns]
    double vd{0.};
    /// Diffusion along E.
    double dSigmaL{0.};
    /// Diffusion transverse to E (radial, phi dir is net 0).
    double dSigmaT{0.};

    double wv{0.};  ///< flux drift velocity [cm/ns]
    double wr{0.};  ///< bulk drift velocity [cm/ns]
    /// Ionization rate from TOF experiment 1/ns -> 1/cm
    double townsendPT{0.};
    /// Attachment rate from TOF experiment 1/ns -> 1/cm
    double attachmentPT{0.};
    /// Space-charge electric field in R direction (can be negative)
    double er{0.};
    /// Space-charge electric field in Z direction (can be negative)
    double ez{0.};

    double time{0.};  ///< Node clock.

    bool anode{false};  ///< init the anode
    /// Gas gap index: -1 if not gas gap; starts with 0, 1, ...
    int gasGapIndex{0};
  };

  struct Point {
    double x{0.};
    double y{0.};
    double z{0.};  ///< coordinates
    double t{0.};  ///< time
  };

  // Prepare grid and place stored electrons from AvalancheMicroscopic import
  void PrepareElectronsFromMicroscopicAvalanche();

  // Assign electron to the closest grid point
  bool SnapTo2dGrid(double x, double y, double z, long n = 1, int gasLayer = 0);

  // Prepare the mesh with the ComponentParallelPlate
  void Prepare2dMesh();

  // Transports the electrons/nodes a timestep
  bool TransportTimeStep();

  // Diffuses the electrons/nodes a timestep
  void DiffuseTimeStep(double dx, long nElectron, double nPosIon,
                       double nNegIon, int iz, int ir, int gasGap);

  // Redistributes the charges
  void DistributeCharges(long nElectron, double nPosIon, double nNegIon, int iz,
                         int ir, double stepZ, double stepR, int gasGap);

  // Get swarm parameters at electric field magnitude
  void GetSwarmParameters(double emag, double &alpha, double &eta,
                          double &vd, double &dSigmaL, double &dSigmaT,
                          double &wv, double &wr, double &alphaPT,
                          double &etaPT, int gasGap);

  // Get from index the gas gap number, else -1
  int GetGasGapNumber(int layerIndex);

  void SetRingSystems();

 private:
  std::string m_className{"AvalancheGridSpaceCharge"};

  bool m_bDebug{false};
  bool m_bDiffusion{false};
  bool m_bStick{true};
  // boolean for AvalancheElectron
  bool m_bDriftAvalanche{false};
  // boolean for ImportElectronsFromAvalancheMicroscopic
  bool m_bImportAvalanche{false};
  bool m_bPreparedImportAvalanche{false};
  long m_lNCrit = {100000000};
  bool m_bSpaceCharge{true};

  float m_fStreamerK{0.95};
  bool m_bStopAtK{false};
  bool m_bFieldK{false};
  long m_lElectronsK{0};

  bool m_bAdaptiveTime{true};
  /// Flag if TOF parameters should be used, else Magboltz
  /// drift and SST spatial coefficients
  bool m_bUseTOF{true};

  bool m_bMC{true};

  int m_iFieldApprox{1};    //< order of approximation in Set(1,2,3,...)

  ComponentParallelPlate *m_pp{nullptr};
  Sensor *m_sensor{nullptr};

  std::vector<double> m_zGrid;  ///< Grid points of z-coordinate.
  int m_zSteps{0};              ///< Number of grid points.
  double m_zStepSize{0.};       ///< Distance between the grid points.
  double m_zInvStep{0.};        ///< Inverse of the grid spacing.

  std::vector<double> m_rGrid;  ///< Grid points of r-coordinate.
  int m_rSteps{0};              ///< Number of grid points
  double m_rStepSize{0.};       ///< Distance between the grid points.
  double m_rInvStep{0.};        ///< Inverse of the grid spacing.

  bool m_isgridset{false};  ///< Keeps track if the grid has been defined.
  long m_nTotElectron{0};   ///< Total amount of electrons at time step.
  long m_nTotPosIons{0};    ///< total amount of charge created

  double m_time{0.};   ///< Clock.
  double m_time0{0.};  ///< Initial time.
  double m_dt{0.};     ///< Time step.

  std::vector<std::vector<int>>
      m_zGasGapBoundaries;  ///< [k] -> {izLeft, ..., izRight}

  std::vector<std::vector<GridNode>> m_grid;  ///< grid with nodes on it
  /// Electrons to transfer onto grid
  std::vector<std::vector<Point>> m_vElectrons;

  std::vector<std::pair<double, long>> m_vNElectronEvolution;
  /// Which layer indices are gas layers
  std::vector<int> m_vIndexGasGaps = {0};

  /// Coordinates of center of electron number.
  /// Required: y in [zmin, zmax]
  std::vector<std::vector<double>> m_vCoNGasLayer;
  /// Example point (y-coord) in each gas gap
  std::vector<double> m_vYPointInGasGap;

  /// Uniform background field in z direction, can be negative.
  std::vector<double> m_ezBkg = {0};
  /// Which gas gaps are saturated if saturation is on
  std::vector<int> m_vSaturatedGaps;

  enum class FieldOption {
    Coulomb,
    Mirror
  };
  FieldOption m_fieldOption{FieldOption::Coulomb};

  /// Vector of ComponentChargedRing objects
  /// We might need multiple ring systems, e.g. one per gas gap.
  std::vector<ComponentChargedRing> m_vRingSystems;
};

}  // namespace Garfield

#endif  // GARFIELD_AVALANCHEGRIDSPACECHARGE_HH
