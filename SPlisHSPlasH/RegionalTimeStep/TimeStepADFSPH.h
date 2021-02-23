#ifndef __TIMESTEPADFSPH_H__
#define __TIMESTEPADFSPH_H__

#include "ParticleGrid.h"
#include "SPlisHSPlasH/TimeStep.h"

namespace SPH
{
    class TimeStepADFSPH : public TimeStep
    {
        public:
            TimeStepADFSPH() = default;
            ~TimeStepADFSPH() = default;

            FORCE_INLINE ParticleGrid& getParticleGrid() { return m_particleGrid; }

        private:
            ParticleGrid m_particleGrid;
    };
}

#endif
