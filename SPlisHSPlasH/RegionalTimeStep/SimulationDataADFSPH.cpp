#include "SimulationDataADFSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

void SimulationDataADFSPH::init(bool isCopy)
{
	m_isCopy = isCopy;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_factor.resize(nModels);
	m_kappa.resize(nModels);
	m_kappaV.resize(nModels);
	m_density_adv.resize(nModels);
	m_oldVelocity.resize(nModels);
	m_averageAcceleration.resize(nModels);
	m_averageAccelerationBorder.resize(nModels);

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_factor[i].resize(fm->numParticles(), 0.0);
		m_kappa[i].resize(fm->numParticles(), 0.0);
		m_kappaV[i].resize(fm->numParticles(), 0.0);
		m_density_adv[i].resize(fm->numParticles(), 0.0);
		m_oldVelocity[i].resize(fm->numParticles());
		m_averageAcceleration[i].resize(fm->numParticles());
		m_averageAccelerationBorder[i].resize(fm->numParticles());
	}
}

void SimulationDataADFSPH::cleanup()
{
	m_factor.clear();
	m_kappa.clear();
	m_kappaV.clear();
	m_density_adv.clear();
	m_oldVelocity.clear();
	m_averageAcceleration.clear();
	m_averageAccelerationBorder.clear();
}

void SimulationDataADFSPH::reset()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numActiveParticles(); j++)
		{
			m_kappa[i][j] = 0.0;
			m_kappaV[i][j] = 0.0;
		}
	}
}

void SimulationDataADFSPH::copyData(const SimulationDataADFSPH* other, unsigned int modelIdx, const unsigned int* particleIndices, const ParticleState* particleStates, int particleCount)
{
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int particleNr = 0; particleNr < particleCount; particleNr++)
		{
			const unsigned int i = particleIndices[particleNr];

			if (particleStates[i] == ParticleState::Active)
			{
				m_factor[modelIdx][i] = other->m_factor[modelIdx][i];
				m_kappa[modelIdx][i] = other->m_kappa[modelIdx][i];
				m_kappaV[modelIdx][i] = other->m_kappaV[modelIdx][i];
				m_density_adv[modelIdx][i] = other->m_density_adv[modelIdx][i];
			}
		}
	}
}

void SimulationDataADFSPH::swapData(SimulationDataADFSPH* other, unsigned int modelIdx, const unsigned int* particleIndices, const ParticleState* particleStates, const unsigned int* isBorder, int particleCount)
{
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int particleNr = 0; particleNr < particleCount; particleNr++)
		{
			const unsigned int i = particleIndices[particleNr];

			if (particleStates[i] == ParticleState::Active && isBorder[i] != UINT32_MAX)
			{
				// Swap data with other data storage
				Real temp = m_factor[modelIdx][i];
				m_factor[modelIdx][i] = other->m_factor[modelIdx][i];
				other->m_factor[modelIdx][i] = temp;

				temp = m_kappa[modelIdx][i];
				m_kappa[modelIdx][i] = other->m_kappa[modelIdx][i];
				other->m_kappa[modelIdx][i] = temp;

				temp = m_kappaV[modelIdx][i];
				m_kappaV[modelIdx][i] = other->m_kappaV[modelIdx][i];
				other->m_kappaV[modelIdx][i] = temp;

				temp = m_density_adv[modelIdx][i];
				m_density_adv[modelIdx][i] = other->m_density_adv[modelIdx][i];
				other->m_density_adv[modelIdx][i] = temp;
			}
		}
	}
}

void SimulationDataADFSPH::clearAvgAccelerations()
{
	Vector3r zeroVec;
	zeroVec.setZero();

	for (std::vector<Vector3r>& accVec : m_averageAcceleration)
		std::fill_n(accVec.data(), accVec.size(), zeroVec);

	for (std::vector<Vector3r>& accVec : m_averageAccelerationBorder)
		std::fill_n(accVec.data(), accVec.size(), zeroVec);
}

void SimulationDataADFSPH::performNeighborhoodSearchSort()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		const unsigned int numPart = fm->numActiveParticles();
		if (numPart != 0)
		{
			auto const& d = sim->getNeighborhoodSearch()->point_set(fm->getPointSetIndex());
			d.sort_field(&m_factor[i][0]);
			d.sort_field(&m_kappa[i][0]);
			d.sort_field(&m_kappaV[i][0]);
			d.sort_field(&m_density_adv[i][0]);
		}
	}
}

void SimulationDataADFSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize kappa values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_kappa[fluidModelIndex][j] = 0.0;
		m_kappaV[fluidModelIndex][j] = 0.0;
	}
}
