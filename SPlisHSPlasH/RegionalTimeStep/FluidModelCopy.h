#ifndef __FluidModelCopy_h__
#define __FluidModelCopy_h__

#include "SPlisHSPlasH/FluidModel.h"

namespace SPH
{
    class NonPressureForceBase;

    /*  ADFSPH requires duplicated particle data to store simulation states at different times. This class keeps such
        duplicated data. */
    class FluidModelCopy : FluidModel
    {
    public:
        FluidModelCopy() = default;
        ~FluidModelCopy() = default;

        // Copy particle data and non-pressure force methods like surface tensions
        void init(const FluidModel* fluidModel);

        void syncParticleData(const FluidModel* fluidModel);

    private:
        void copyNonPressureForceMethods(const FluidModel* fluidModel);
        template <typename T>
        void copyNonPressureForceMethod(int methodTypeID, T** methodBase, unsigned int& currentMethodNr, unsigned int newMethodNr);
        void copySurfaceTensionMethod(const FluidModel* fluidModel);
    };

    template <typename T>
    void FluidModelCopy::copyNonPressureForceMethod(int methodTypeID, T** methodBase, unsigned int& currentMethodNr, unsigned int newMethodNr)
    {
        Simulation* sim = Simulation::getCurrent();
        const std::vector<Simulation::NonPressureForceMethod>& stMethods = sim->getSurfaceTensionMethods();

        if (newMethodNr >= stMethods.size())
            newMethodNr = 0;

        if (newMethodNr == m_surfaceTensionMethod)
            return;

        delete *methodBase;
        *methodBase = nullptr;

        currentMethodNr = newMethodNr;
        for (unsigned int i = 0; i < stMethods.size(); i++)
        {
            if (stMethods[i].m_id == methodTypeID)
                *methodBase = static_cast<T*>(stMethods[i].m_creator(this));
        }
    }
}

#endif
