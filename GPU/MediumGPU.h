#ifndef G_MEDIUMGPU_H
#define G_MEDIUMGPU_H

#ifndef __GPUCOMPILE__
#error GPU HEADER INCLUDED WITHOUT SETTING __GPUCOMPILE__
#endif

#include "Garfield/HelperMacros.hh"
#include "GPUInterface.hh"
#include "Garfield/MagboltzInterface.hh"
#include "Garfield/GarfieldConstants.hh"

namespace Garfield
{

class MediumGPU
{
public:
    /// Constructor
  MediumGPU() = default;
  /// Destructor
  ~MediumGPU() {};
  /// Return the id number of the class instance.
  __DEVICE__ int GetId() const { return m_id; }
  /// Is charge carrier transport enabled in this medium?
  __DEVICE__ bool IsDriftable() const { return m_driftable; }
  /// Does the medium have electron scattering rates?
  __DEVICE__ bool IsMicroscopic() const { return m_microscopic; }
  __device__ cuda_t GetElectronCollisionRate(const cuda_t e, const int band);

  __device__ bool ElectronCollision(const cuda_t e, int& type, int& level,
                                    cuda_t& e1, cuda_t& dx, cuda_t& dy,
                                    cuda_t& dz, Particle* secondaries_type,
                                    cuda_t* secondaries_energy,
                                    int& num_secondaries, int& ndxc, int& band);
    // Id number
  int m_id;

  // Transport flags
  bool m_driftable = false;
  bool m_microscopic = false;
  bool m_ionisable = false;
  #include "Garfield/MediumMagboltz.hh"

  friend class MediumGas;
  friend class MediumMagboltz;

  // enum to mimic polymorphism
  enum class MediumType { Medium = 0, MediumGas, MediumMagboltz };

  MediumType m_MediumType{MediumType::Medium};
private:
};

}
#endif
