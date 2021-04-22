#ifndef __TIMESTEPADFSPH_H__
#define __TIMESTEPADFSPH_H__

#include "FluidModelCopy.h"
#include "ParticleGrid.h"
#include "SimulationDataADFSPH.h"
#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/RegionalTimeStep/ParticleGrid.h"

#define USE_WARMSTART
#define USE_WARMSTART_V

namespace SPH
{
	class TimeStepADFSPH : public TimeStep
	{
	protected:
		SimulationDataADFSPH m_simulationData;
		SimulationDataADFSPH m_simulationDataCopy;
		ParticleGrid m_particleGrid;
		unsigned int m_counter;
		const Real m_eps = static_cast<Real>(1.0e-5);
		bool m_enableDivergenceSolver;
		unsigned int m_iterationsV;
		Real m_maxErrorV;
		unsigned int m_maxIterationsV;
		/*	m_subStepNr helps keep track of which regional levels need to be calculated in the current step.
			Has a value within the interval [0, TIMESTEP_MULTIPLIER^(LEVEL_COUNT - 1)). */
		unsigned int m_subStepNr;
		unsigned int m_lastCalculatedLevel;
		unsigned int m_highestLevelToStep; // The highest region level that is or has been stepped in the current step

		std::vector<FluidModelCopy*> m_fluidModelCopies;
		std::vector<FluidModel*> m_originalFluidModel;

		void checkVelocities();
		void computeDFSPHFactor(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void divergenceSolve();
		void divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void computeDensityAdv(unsigned int fluidModelIndex, unsigned int particleNr, unsigned int i, Real h, Real density0);
		void computeDensityChange(unsigned int fluidModelIndex, unsigned int particleNr, unsigned int i, Real h);

#ifdef USE_WARMSTART_V
		void warmstartDivergenceSolve(const unsigned int fluidModelIndex);
#endif
#ifdef USE_WARMSTART
		void warmstartPressureSolve(const unsigned int fluidModelIndex);
#endif

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();
		void calculateLevel(unsigned int level, Real dt);
		void interpolateBorderParticles(unsigned int level);
		void correctAccelerations(unsigned int level);

		void setActiveParticles(unsigned int level);

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		virtual void initParameters();

		void debugParticle(unsigned int fluidModelIndex, unsigned int i);

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
