#ifndef G_RANDOM_ENGINE_PRECALCGPU_H
#define G_RANDOM_ENGINE_PRECALCGPU_H

#include "GPUInterface.hh"

namespace Garfield {

class RandomEnginePreCalcGPU {
 public:
  /// Constructor
  RandomEnginePreCalcGPU();
  /// Destructor
  ~RandomEnginePreCalcGPU();
  /// Call the random number generator.
  __device__ GPUFLOAT Draw();
  void setRandomEngineOnDevice();

// private:
  GPUFLOAT **m_draws;
  unsigned int *m_currIdx;
  unsigned int m_maxSeed{0};
  unsigned int m_maxIdx{0};

  friend class RandomEnginePreCalc;
};

}

#endif
