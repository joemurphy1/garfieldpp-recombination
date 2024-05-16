#include "RandomEnginePreCalcGPU.h"
#include "Garfield/RandomEnginePreCalc.hh"
#include "GPUFunctions.h"
#include <iostream>

namespace Garfield {
#ifdef USEPRECALCRNG
    __device__ RandomEnginePreCalcGPU *randomEngineGPU;

    static constexpr double Pi = 3.1415926535897932384626433832795;
    static constexpr double TwoPi = 2. * Pi;

    __global__ void setRandomEngine_d(RandomEnginePreCalcGPU *engine)
    {
        randomEngineGPU = engine;
    }

    __device__ GPUFLOAT RndmUniformGPU()
    {
        return randomEngineGPU->Draw();
    }

    __device__ GPUFLOAT RndmUniformPosGPU() {
        GPUFLOAT r = RndmUniformGPU();
        while (r <= 0.) r = RndmUniformGPU();
        return r;
    }

    __device__ void RndmDirectionGPU(GPUFLOAT& dx, GPUFLOAT& dy, GPUFLOAT& dz,
        const GPUFLOAT length) {
        const GPUFLOAT phi = TwoPi * RndmUniformGPU();
        const GPUFLOAT ctheta = 2 * RndmUniformGPU() - 1.;
        const GPUFLOAT stheta = sqrt(1. - ctheta * ctheta);
        dx = length * cos(phi) * stheta;
        dy = length * sin(phi) * stheta;
        dz = length * ctheta;
    }

    __device__ void ResetSeed() {
        randomEngineGPU->m_currIdx[(threadIdx.x + blockIdx.x * blockDim.x)] = 0;
    }

#endif

    RandomEnginePreCalcGPU::RandomEnginePreCalcGPU(){}

    RandomEnginePreCalcGPU::~RandomEnginePreCalcGPU()
    {
    }

    void RandomEnginePreCalcGPU::setRandomEngineOnDevice()
    {
#ifdef USEPRECALCRNG
        setRandomEngine_d<<<1,1>>>(this);
#endif
    }

    __device__ GPUFLOAT RandomEnginePreCalcGPU::Draw()
    {
        unsigned int seed = (threadIdx.x + blockIdx.x * blockDim.x);
        double draw{m_draws[seed][m_currIdx[seed]++]};
#ifndef GPUOPTIMISE
        if (m_currIdx[seed] >= MAXRNGDRAWS)
        {
            m_currIdx[seed] = 0;
        }
#endif
        return draw;

    }

    double RandomEnginePreCalc::CreateGPUTransferObject(RandomEnginePreCalcGPU *&rng_gpu)
    {
        checkCudaErrors(cudaMallocManaged(&rng_gpu, sizeof(RandomEnginePreCalcGPU)));

        rng_gpu->m_maxSeed = m_maxSeed;
        rng_gpu->m_maxIdx = m_maxIdx;

        checkCudaErrors( cudaMalloc(&rng_gpu->m_draws, m_maxSeed*sizeof(GPUFLOAT*)));
        checkCudaErrors( cudaMalloc(&rng_gpu->m_currIdx, m_maxSeed * sizeof(int)));

        // assign memory for temporary pointer store
        GPUFLOAT **h_rng_array_ptrs = new GPUFLOAT*[m_maxSeed];
        for (unsigned int i = 0; i < m_maxSeed; i++)
        {
            checkCudaErrors( cudaMalloc(&(h_rng_array_ptrs[i]), m_maxIdx*sizeof(GPUFLOAT)));
        }

        // copy from host to device memory
        checkCudaErrors( cudaMemcpy(rng_gpu->m_draws, h_rng_array_ptrs, m_maxSeed*sizeof(GPUFLOAT*), cudaMemcpyHostToDevice) );

        for (int i = 0; i < m_maxSeed; i++)
        {
            checkCudaErrors( cudaMemcpy(h_rng_array_ptrs[i], &(m_draws[i][0]), m_maxIdx*sizeof(GPUFLOAT), cudaMemcpyHostToDevice));
        }

        return sizeof(RandomEnginePreCalcGPU) + m_maxSeed*sizeof(GPUFLOAT*) +
            m_maxSeed * sizeof(int) + m_maxSeed * m_maxIdx * sizeof(GPUFLOAT);

    }
}
