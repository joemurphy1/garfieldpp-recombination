#include "Garfield/RandomEnginePreCalc.hh"
#include <iostream>

#include "GPUInterface.hh"

namespace Garfield {

#ifdef USEPRECALCRNG
// only create this if we're using it as it takes a while to fill
RandomEnginePreCalc randomEngine;
#endif

RandomEnginePreCalc::RandomEnginePreCalc() : RandomEngine() {
    std::cout << "RandomEnginePreCalc:  Populating RNG..." << std::endl;
    m_draws.resize(m_maxSeed);
    m_currIdx.resize(m_maxSeed);

    for (unsigned int i = 0; i < m_maxSeed; i++)
    {
        m_currIdx[i] = 0;
        m_rng.SetSeed(i + 1);
        m_draws[i].resize(m_maxIdx);

        for(unsigned int j = 0; j < m_maxIdx; j++)
        {
            m_draws[i][j] = m_rng.Rndm();
        }
    }

    std::cout << "RandomEnginePreCalc:  RNG Populated" << std::endl;
}

RandomEnginePreCalc::~RandomEnginePreCalc() 
{
}

double RandomEnginePreCalc::Draw() 
{
    double draw{m_draws[m_seed][m_currIdx[m_seed]++]};
    if (m_currIdx[m_seed] >= m_draws[m_seed].size())
    {
        std::cout << "WARNING: Wrapping draw from seed " << m_seed << ". Reached " << m_draws[m_seed].size() << std::endl;
        m_currIdx[m_seed] = 0;
    }
    return draw;
}

void RandomEnginePreCalc::Seed(const unsigned int s) {

    m_seed = s % m_maxSeed;
    m_currIdx[m_seed] = 0;

    if (s >= m_maxSeed)
    {
        std::cout << "WARNING: Wrapping seed. Requested " << m_maxSeed << std::endl;
    }
}

void RandomEnginePreCalc::Print() {
  std::cout << "RandomEngineRoot::Print:\n"
            << "    Generator type: TRandom3\n"
            << "    Seed: " << 1 << "\n";
}

}
