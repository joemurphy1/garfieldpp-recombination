#ifndef G_RANDOMGPU_H
#define G_RANDOMGPU_H

#include "GPUInterface.hh"
#include "RandomEnginePreCalcGPU.h"

namespace Garfield {

  __device__ void ResetSeed();
  __device__ GPUFLOAT RndmUniformGPU();
  __device__ GPUFLOAT RndmUniformPosGPU();
  __device__ void RndmDirectionGPU(GPUFLOAT& dx, GPUFLOAT& dy, GPUFLOAT& dz, const GPUFLOAT length = 1.);

}

#endif