#ifndef GARFIELD_AVALANCHEGRIDSPACECHARGE_HH
#define GARFIELD_AVALANCHEGRIDSPACECHARGE_HH

#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <numeric>

#include "AvalancheMicroscopic.hh"
#include "ComponentParallelPlate.hh"
#include "ComponentConstant.hh"
#include "GarfieldConstants.hh"
#include "Sensor.hh"
#include "Random.hh"


// Get size of avalanche when going from x to x+dx in Monte Carlo fashion
void GetAvalancheSizeFromStep(double dx, const long nElectronIn, const double alpha, const double eta,
                              long &nElectronOut, double &nPosIonOut, double &nNegIonOut);

// Get mean size of avalanche when going from x to x+dx
void GetMeanAvalancheSizeFromStep(double dx, const long nElectronIn, const double alpha, const double eta,
                                  long &nElectronOut, double &nPosIonOut, double &nNegIonOut);


namespace Garfield {
  /// Propagates avalanches with the 2d (axi-symmetric) space-charge routine from Lippmann, Riegler (2004)
  /// in uniform background fields. Different options to calculate space-charge-fields can be chosen.
  class AvalancheGridSpaceCharge {
  public:
    /// Constructor
    AvalancheGridSpaceCharge();

    /// Destructor
    ~AvalancheGridSpaceCharge() = default;

    /// Reset grid i.e. ImportElectrons can be called again.
    void Reset();

    /// Enable/disable debugging ( = log messages) (default off)
    void EnableDebugging(const bool option = true) { m_bDebug = option; }

    /// Enable sticky anode (default on)
    void EnableStickyAnode(const bool option = true) { m_bStick = option; }

    /// Enable to use TOF swarm parameters (default on)
    void EnableTOF(const bool option = true) { m_bUseTOF = option; }

    /// Enable diffusion (default off)
    void EnableDiffusion(const bool option = true) { m_bDiffusion = option; }

    /// Enable space charge calculations (default on)
    void EnableSpaceChargeEffect(const bool option = true) { m_bSpaceCharge = option; }

    /// Enable adaptive time stepping (default on)
    void EnableAdaptiveTimeStepping(const bool option = true) { m_bAdaptiveTime = option; }

    /// Enable Monte Carlo (gain/diffusion) up to a certain total number of electrons e.g. (1e5) (default on)
    ///   if disabled: Mean values are considered instead of the statistical processes
    void EnableMC(const bool option = true) { m_bMC = option; }

    /// Set electron saturation value applied for each gas gap if space charge effect is turned off (default 1e8)
    void SetNCrit(const long NCrit = 1e8) { m_lNCrit = NCrit; }

    /// Sets the different options to calculate the space charge field
    ///   coulomb: free field approximation
    ///   (relaxation: fdm-relaxation methode (with a condition to stop))
    ///   mirror: symmetric 3 layer single gap rpc with metal - resistive layer - gas gap - r. l. - m.
    void SetFieldCalculation(const std::string &option = "coulomb", const int nof_approx = 1.) {
      m_sFieldOption = std::move(option);
      m_iFieldApprox = std::move(nof_approx);
    }

    /// Set the streamer-inception criterion constant K in the interval (0, \infty) s.t. 1 = 100%
    void SetK(float option = 0.95) { m_fStreamerK = option; }

    /// Stop the avalanche if K % field is reached
    void SetStopAtK(bool option = true) { m_bStopAtK = option; }

    /// Set the sensor (+ determines if base Cmp is CmpParallelPlate (MRPCS)).
    void SetSensor(Sensor *sensor) {
      // set Sensor
      m_sensor = sensor;
      // determine if any component is CmpParallelPlate (if not it will stay nullptr)
      size_t nofCmp = m_sensor->GetNumberOfComponents();
      for (int i=0; i < nofCmp; i++) {
        if (!m_ParallelPlate) {
          m_ParallelPlate = dynamic_cast<ComponentParallelPlate *>(m_sensor->GetComponent(i));
        }
      }
    }

    /**
     *
     * @param zmin coordinate along direction of background field
     * @param zmax should more or less be at the end of the gap
     * @param zsteps number of steps
     * @param rmax coordinate perpendicular to background field
     * @param rsteps number of steps
     */
    void Set2dGrid(double zmin, double zmax, int zsteps,
                   double rmax, int rsteps);

    /// Import electron (no ions) data from AvalancheMicroscopic class to axi-symmetric grid
    void ImportElectronsFromAvalancheMicroscopic(AvalancheMicroscopic *avmc);

    /// Set n electrons onto the grid
    void AvalancheElectron(double x, double y, double z,
                           double t = 0, int n = 1);

    /// (After calling AvalancheElectron) add more electrons on the same transversal line (y-freedom)
    void AddExtraAvalancheElectron(double y, int n = 1);

    /// Starts the simulation with the imported electrons for a time step dt
    /// (dt = -1 until there are no electrons left in the gap).
    void StartGridAvalanche(double dtime = -1);

    /// Returns the total positive charge in the gap's
    long GetAvalancheSize() { return m_AvGrid.nTotPosIons; }

    /// Return current mean distance of the electrons on the grid
    double GetMeanDistance();

    /// Returns if 100 * K % background charge has been reached
    [[nodiscard]] long ReachedKPercent() const {
      if (m_bFieldK) return m_lElectronsK;
      else return -1;
    }

    /// Returns the total electron number evolution
    [[nodiscard]] const std::vector<std::pair<double, long>> GetElectronEvolution() const {
      return m_vNElectronEvolution;
    }

    /// Export the current grid to a txt file (electron, ions numbers and field magnitude)
    ///  Take care: only where electrons are located is the space.charge field evaluated!
    void ExportGrid(const std::string &filename);

  private:
    struct GridNode {
      long nElectron = 0; ///< electrons on node
      double nPosIon = 0; ///< pos ion on node (smeared values allowed)
      double nNegIon = 0; ///< neg ion on node (smeared values allowed)
      // holder memories for stepping in time:
      long nElectronHolder = 0; ///< at t+dt
      double nPosIonHolder = 0; ///< at t+dt
      double nNegIonHolder = 0; ///< at t+dt

      double townsend = 0; ///< townsend at this node 1/cm
      double attachment = 0; ///< attachment at this node 1/cm
      double velocity = 0; ///< magnitude of velocity of the node (not negative) cm/ns
      double dSigmaL = 0; ///< diffusion along E
      double dSigmaT = 0; ///< diffusion transverse to E (radial, phi dir is netto 0)
      double Wv = 0; ///< flux drift cm/ns
      double Wr = 0; ///< bulk drift cm/ns
      double townsendPT = 0; ///< ionization rate from TOF experiment 1/ns -> 1/cm
      double attachmentPT = 0; ///< attachment rate from TOF experiment 1/ns -> 1/cm

      double eFieldR = 0; ///< space-charge electric field in R direction (can be negative)
      double eFieldZ = 0; ///< space-charge electric field in Z direction (can be negative)

      double time = 0.;  ///< Node clock.

      bool anode = false; ///< init the anode

      int layerIndex = 0; ///< layerIndex in ParallelPlate convention != gas gap index
      int gasGapIndex = 0; ///< gas gap index: -1 if not gas gap; starts with 0, 1, ...
      bool isGasGap = true;
    };

    struct Grid {
      std::vector<double> zGrid;  ///< Grid points of z-coordinate.
      int zSteps = 0;             ///< Amount of grid points.
      double zStepSize = 0.;      ///< Distance between the grid points of z-coordinate.

      std::vector<std::vector<int>> zGasGapBoundaries; ///< [k] -> {izLeft, ..., izRight}

      std::vector<double> rGrid;  ///< Grid points of x-coordinate.
      double rStepSize = 0.;      ///< Amount of grid points.
      int rSteps = 0.;            ///< Distance between the grid points of x-coordinate.

      bool isgridset = false;     ///< Keeps track if the grid has been defined.
      long nTotElectron = 0;      ///< Total amount of electrons at time step.
      long nTotPosIons = 0;       ///< total amount of charge created

      double time = 0.;           ///< Grid clock.
      double time0 = 0.;          ///< initial time
      double dt = 0.;             ///< time step.

      bool run = true;            ///< Tracking if the charges are still in the drift gap.
    };

    struct Point {
      double x, y, z; ///< coordinates
      double t; ///< time

      int gasLayerIndex;
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
    void DiffuseTimeStep(double dx, long nElectron, double nPosIon, double nNegIon,
                         int iz, int ir, int gasGap);

    // Redistributes the charges
    void DistributeCharges(long nElectron, double nPosIon, double nNegIon,
                           int iz, int ir, double stepZ, double stepR, int gasGap);

    // Calculate the field from all the contributions to the bin of interest. May need much more functionalities/tables.
    void GetLocalField(int iz, int ir, double &eFieldZ, double &eFieldR, const std::string &fieldOption,
                       int gasGap);

    // Calculate the field of charged ring in vacuum using coulomb potentials and indices
    void GetFreeChargedRing(int iz, int ir, int fz, int fr, double &eFieldZ, double &eFieldR);

    // Calculate the field of charged ring in vacuum using coulomb potentials and coordinates
    void GetFreeChargedRing(double zi, double ri, double zf, double rf, double &eFieldZ, double &eFieldR);

    // Get field at (zi, ri) from N charges at (zf, rf) either as a ring or a coulomb ball (rf = 0)
    // if i and f are too close it is considered as self interaction and not included
    bool AddFieldFromChargeAt(int iz, int ir, int fz, int fr, double N, double &eFieldZ, double &eFieldR);

    // Get field at (zi, ri) from N charges at (zf, rf) either as a ring or a coulomb ball (rf = 0)
    // if i and f are too close it is considered as self interaction and not included
    bool AddFieldFromChargeAt(int iz, int ir, double zf, double rf, double N, double &eFieldZ, double &eFieldR);

    // Get swarm parameters at electric field magnitude
    void GetSwarmParameters(double MagEField, double &alpha, double &eta, double &drift,
                            double &dSigmaL, double &dSigmaT, double &wv, double &wr,
                            double &alphaPT, double &etaPT, int gasGap);

    // Change from 2dGrid to Global coordinates
    void GetGlobalCoordinates(double r, double z, double phi, double &xg, double &yg, double &zg, int gasGap);

    // Import elliptic integral values
    void ImportEllipticIntegralValues(const std::string &filename);

    // Gets elliptic integrals via list
    void GetEllipticIntegrals(double x, double &K, double &E);

    // Get from index the gas gap number, else -1
    int GetGasGapNumber(int layerIndex) {
      auto it = std::find(m_vIndexGasGaps.begin(), m_vIndexGasGaps.end(), layerIndex);
      return (it != m_vIndexGasGaps.end()) ? std::distance(m_vIndexGasGaps.begin(), it) : -1;
    }

  private:
    std::string m_className = "AvalancheGridSpaceCharge";

    bool m_bDebug = false;
    bool m_bDiffusion = false;
    bool m_bStick = true;
    // boolean for AvalancheElectron
    bool m_bDriftAvalanche = false;
    // boolean for ImportElectronsFromAvalancheMicroscopic
    bool m_bImportAvalanche = false;
    bool m_bPreparedImportAvalanche = false;
    long m_lNCrit = 1e8;
    bool m_bSpaceCharge = true;

    float m_fStreamerK = 0.95;
    bool m_bStopAtK = false;
    bool m_bFieldK = false;
    long m_lElectronsK;

    bool m_bAdaptiveTime = true;
    bool m_bImportElliptic = false;
    bool m_bUseTOF = true; //< if TOF parameters should be used else Magboltz drift and SST spatial coefficients
    bool m_bWrAvailable = true; //< if bulk drift velocity is available to the simulation
    bool m_bRatesAvailable = true; //< if temporal rates are available to the simulation
    bool m_bMC = true;

    int m_iFieldApprox = 1; //< order of approximation in Set(1,2,3,...)
    double m_dMinGroups = 50; //< same values as lippmann

    ComponentParallelPlate *m_ParallelPlate = nullptr;
    Sensor *m_sensor = nullptr;
    Grid m_AvGrid;
    std::vector<std::vector<GridNode>> m_GridMesh; ///< grid with nodes on it
    std::vector<std::vector<Point>> m_vElectrons; ///< reminder of electrons to transfer onto grid
    std::vector<std::pair<double, long>> m_vNElectronEvolution;
    std::vector<long> m_vGroupSizes = {1500, 800, 400, 200, 100, 50, 20, 10, 5, 2}; ///< same values as Lippmann & Riegler
    std::vector<int> m_vIndexGasGaps = {0}; ///< which layer indices are gas layers
    std::vector<std::vector<double>> m_vCoNGasLayer{}; ///< coordinates of center of electron number. Required: y in [zmin, zmax]
    std::vector<double> m_vYPointInGasGap{}; ///< example point (y-coord) in each gas gap
    std::vector<double> m_vEFieldZBackgroundGasLayer = {
            0}; ///< uniform background field in z direction, can be negative.
    std::vector<int> m_vSaturatedGaps{}; ///< which gas gaps are saturated if saturation is on
    std::string m_sFieldOption = "coulomb";
    std::vector<double> m_vXElliptic;
    std::vector<double> m_vKElliptic;
    std::vector<double> m_vEElliptic;
  };

} // namespace Garfield

#endif //GARFIELD_AVALANCHEGRIDSPACECHARGE_HH
