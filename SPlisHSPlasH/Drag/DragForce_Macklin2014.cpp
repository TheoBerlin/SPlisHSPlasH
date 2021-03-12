#include "DragForce_Macklin2014.h"
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;

DragForce_Macklin2014::DragForce_Macklin2014(FluidModel *model) :
	DragBase(model)
{
}

DragForce_Macklin2014::~DragForce_Macklin2014(void)
{
}

void DragForce_Macklin2014::step()
{
	const Real density0 = m_model->getDensity0();

	const int numParticles = (int)m_model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		const unsigned int* particleIndices = m_model->getParticleIndices();

		#pragma omp for schedule(static)
		for (int particleNr = 0; particleNr < numParticles; particleNr++)
		{
			const unsigned int i = particleIndices[particleNr];
			Vector3r &ai = m_model->getAcceleration(i);
			const Vector3r &vi = m_model->getVelocity(i);
			ai -= m_dragCoefficient * static_cast<Real>(1.0) / m_model->getMass(i) * vi * (1.0 - m_model->getDensity(i) / density0);
		}
	}
}


void DragForce_Macklin2014::reset()
{
}

