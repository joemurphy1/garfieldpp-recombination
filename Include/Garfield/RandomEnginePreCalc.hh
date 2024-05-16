#ifndef G_RANDOM_ENGINE_PRECALC_H
#define G_RANDOM_ENGINE_PRECALC_H

#include <TRandom3.h>

#include "RandomEngine.hh"
#include "GPUInterface.hh"

namespace Garfield {

class RandomEnginePreCalcGPU;

class RandomEnginePreCalc : public RandomEngine {
 public:
  /// Constructor
  RandomEnginePreCalc();
  /// Destructor
  ~RandomEnginePreCalc();
  /// Call the random number generator.
  double Draw() override;
  /// Initialise the random number generator.
  void Seed(const unsigned int s) override;
  /// Print information about the generator used and the seed. 
  void Print() override;
  /// Retrieve the seed that was used
  unsigned int GetSeed() override {return 0;}

  /// Create and initialise GPU Transfer class
  double CreateGPUTransferObject(RandomEnginePreCalcGPU *&rng_gpu);

 private:
  std::vector< std::vector<GPUFLOAT> > m_draws;
  std::vector<unsigned int> m_currIdx;
  unsigned int m_seed{0};
  unsigned int m_maxSeed{MAXSTACKSIZE};
  unsigned int m_maxIdx{MAXRNGDRAWS};

  TRandom3 m_rng;
};
}

#endif
