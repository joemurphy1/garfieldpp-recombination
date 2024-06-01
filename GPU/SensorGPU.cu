#define __GPUCOMPILE__
#include "Sensor.cc"
#undef __GPUCOMPILE__

#include "Garfield/Sensor.hh"

namespace Garfield {

  double Sensor::CreateGPUTransferObject(SensorGPU *&sensor_gpu)
  {
      // create main sensor GPU class
      checkCudaErrors(cudaMallocManaged(&sensor_gpu, sizeof(SensorGPU)));
      double alloc{sizeof(SensorGPU)};

      // transfer sizes
      sensor_gpu->m_xMinUser = m_xMinUser;
      sensor_gpu->m_yMinUser = m_yMinUser;
      sensor_gpu->m_zMinUser = m_zMinUser;
      sensor_gpu->m_xMaxUser = m_xMaxUser;
      sensor_gpu->m_yMaxUser = m_yMaxUser;
      sensor_gpu->m_zMaxUser = m_zMaxUser;

      // TN: b332e924 introduced a change in m_components. It now is a vector of
      // std::pair<Component*, bool>, which allows the option to disable
      // components. This is NOT supported on the GPU for now, so I make a new
      // vector just of Component* and cross my fingers this will still work
      std::vector<Component*> components{};
      for (const auto& cmp: m_components) {
        components.push_back(std::get<0>(cmp));
      }

      // create arrays
      alloc += CreateGPUObjectArrayFromVector<Component*, ComponentGPU**>(components,
        sensor_gpu->m_numComponents, sensor_gpu->m_components);

      return alloc;
  }
}