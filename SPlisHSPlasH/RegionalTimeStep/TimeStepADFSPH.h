#ifndef __TIMESTEPADFSPH_H__
#define __TIMESTEPADFSPH_H__

#include "FluidModelCopy.h"
#include "ParticleGrid.h"
#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/DFSPH/SimulationDataDFSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/RegionalTimeStep/ParticleGrid.h"

#define USE_WARMSTART
#define USE_WARMSTART_V

namespace SPH
{
	class TimeStepADFSPH : public TimeStep
	{
	protected:
		SimulationDataDFSPH m_simulationData;
		ParticleGrid m_particleGrid;
		unsigned int m_counter;
		const Real m_eps = static_cast<Real>(1.0e-5);
		bool m_enableDivergenceSolver;
		unsigned int m_iterationsV;
		Real m_maxErrorV;
		unsigned int m_maxIterationsV;

		std::vector<FluidModelCopy*> m_fluidModelCopies;

		void checkVelocities();
		void computeDFSPHFactor(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void divergenceSolve();
		void divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int index, const Real h, const Real density0);
		void computeDensityChange(const unsigned int fluidModelIndex, const unsigned int index, const Real h);

#ifdef USE_WARMSTART_V
		void warmstartDivergenceSolve(const unsigned int fluidModelIndex);
#endif
#ifdef USE_WARMSTART
		void warmstartPressureSolve(const unsigned int fluidModelIndex);
#endif

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();
		void stepLevel(unsigned int level, Real dt);

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		virtual void initParameters();

	public:
		static int SOLVER_ITERATIONS_V;
		static int MAX_ITERATIONS_V;
		static int MAX_ERROR_V;
		static int USE_DIVERGENCE_SOLVER;
		static int RENDER_REGION_COLORS;

		TimeStepADFSPH();
		virtual ~TimeStepADFSPH(void);

		virtual void step();
		virtual void reset();

		virtual void resize();

        FORCE_INLINE ParticleGrid& getParticleGrid() { return m_particleGrid; }
	};
}

#endif
