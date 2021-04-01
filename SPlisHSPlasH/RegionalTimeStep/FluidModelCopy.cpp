#include "FluidModelCopy.h"

#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Viscosity/ViscosityBase.h"
#include "SPlisHSPlasH/SurfaceTension/SurfaceTensionBase.h"
#include "SPlisHSPlasH/Vorticity/VorticityBase.h"
#include "SPlisHSPlasH/Drag/DragBase.h"
#include "SPlisHSPlasH/Elasticity/ElasticityBase.h"

using namespace SPH;

void FluidModelCopy::init(const FluidModel* fluidModel)
{
    m_id = fluidModel->m_id + " copy";

    const unsigned int numParticles = fluidModel->m_x.size();

    releaseFluidParticles();
    resizeFluidParticles(numParticles);
    copyParticleData(fluidModel);
    copyNonPressureForceMethods(fluidModel);
}

void FluidModelCopy::copyNonPressureForceMethods(const FluidModel* fluidModel)
{
    copyNonPressureForceMethod(fluidModel->getValue<int>(FluidModel::SURFACE_TENSION_METHOD), &m_surfaceTension,
        m_surfaceTensionMethod, fluidModel->m_surfaceTensionMethod);

    copyNonPressureForceMethod(fluidModel->getValue<int>(FluidModel::VISCOSITY_METHOD), &m_viscosity,
        m_viscosityMethod, fluidModel->m_viscosityMethod);

    copyNonPressureForceMethod(fluidModel->getValue<int>(FluidModel::VORTICITY_METHOD), &m_vorticity,
        m_vorticityMethod, fluidModel->m_vorticityMethod);

    copyNonPressureForceMethod(fluidModel->getValue<int>(FluidModel::DRAG_METHOD), &m_drag,
        m_dragMethod, fluidModel->m_dragMethod);

    copyNonPressureForceMethod(fluidModel->getValue<int>(FluidModel::ELASTICITY_METHOD), &m_elasticity,
        m_elasticityMethod, fluidModel->m_elasticityMethod);
}
