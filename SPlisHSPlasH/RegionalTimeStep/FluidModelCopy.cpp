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

    releaseFluidParticles();
    resizeFluidParticles(fluidModel->m_x.size());
    syncParticleData(fluidModel);
    copyNonPressureForceMethods(fluidModel);
}

void FluidModelCopy::syncParticleData(const FluidModel* fluidModel)
{
    m_masses = fluidModel->m_masses;
    m_a = fluidModel->m_a;
    m_v0 = fluidModel->m_v0;
    m_x0 = fluidModel->m_x0;
    m_x = fluidModel->m_x;
    m_x = fluidModel->m_x;
    m_v = fluidModel->m_v;
    m_density = fluidModel->m_density;
    m_particleId = fluidModel->m_particleId;
    m_particleState = fluidModel->m_particleState;
    m_V = fluidModel->m_V;

#ifdef USE_PERFORMANCE_OPTIMIZATION
    m_precomp_V_gradW = fluidModel->m_precomp_V_gradW;
    m_precompIndices = fluidModel->m_precompIndices;
    m_precompIndicesSamePhase = fluidModel->m_precompIndicesSamePhase;
#endif

    m_density0 = fluidModel->m_density0;
    m_pointSetIndex = fluidModel->m_pointSetIndex;
    m_numActiveParticles = fluidModel->m_numActiveParticles;
    m_numActiveParticles0 = fluidModel->m_numActiveParticles0;
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
