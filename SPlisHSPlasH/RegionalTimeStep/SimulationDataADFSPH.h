#ifndef __SimulationDataADFSPH_h__
#define __SimulationDataADFSPH_h__

#include "SPlisHSPlasH/DFSPH/SimulationDataDFSPH.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"
#include "ParticleGrid.h"

namespace SPH
{
	/** \brief Simulation data which is required by the method Divergence-free Smoothed Particle Hydrodynamics introduced
	* by Bender and Koschier [BK15,BK17].
	*
	* References:
	* - [BK15] Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In ACM SIGGRAPH / Eurographics Symposium on Computer Animation, SCA '15, 147-155. New York, NY, USA, 2015. ACM. URL: http://doi.acm.org/10.1145/2786784.2786796
	* - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
	*/
	class SimulationDataADFSPH
	{
		public:
			SimulationDataADFSPH() = default;
			virtual ~SimulationDataADFSPH() = default;

		protected:
			bool m_isCopy;

			/** \brief factor \f$\alpha_i\f$ */
			std::vector<std::vector<Real>> m_factor;
			/** \brief stores \f$\kappa\f$ value of last time step for a warm start of the pressure solver */
			std::vector<std::vector<Real>> m_kappa;
			/** \brief stores \f$\kappa^v\f$ value of last time step for a warm start of the divergence solver */
			std::vector<std::vector<Real>> m_kappaV;
			/** \brief advected density */
			std::vector<std::vector<Real>> m_density_adv;
			// Velocity of each particle before pressure solving. Used to calculate acceleration.
			std::vector<std::vector<Vector3r>> m_oldVelocity;
			std::vector<std::vector<Vector3r>> m_averageAcceleration;
			std::vector<std::vector<Vector3r>> m_averageAccelerationBorder;

		public:

			/** Initialize the arrays containing the particle data.
			*/
			virtual void init(bool isCopy);

			/** Release the arrays containing the particle data.
			*/
			virtual void cleanup();

			/** Reset the particle data.
			*/
			virtual void reset();

			void copyData(const SimulationDataADFSPH* other, unsigned int modelIdx, const unsigned int* particleIndices, const ParticleState* particleStates, int particleCount);

			// Swaps data with other data storage if a particle is part of a border
			void swapData(SimulationDataADFSPH* other, unsigned int modelIdx, const unsigned int* particleIndices, const ParticleState* particleStates, const unsigned int* isBorder, int particleCount);

			void clearAvgAccelerations();

			/** Important: First call m_model->performNeighborhoodSearchSort()
			 * to call the z_sort of the neighborhood search.
			 */
			void performNeighborhoodSearchSort();
			void emittedParticles(FluidModel *model, const unsigned int startIndex);

			FORCE_INLINE const Real getFactor(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_factor[fluidIndex][i];
			}

			FORCE_INLINE Real& getFactor(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_factor[fluidIndex][i];
			}

			FORCE_INLINE void setFactor(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_factor[fluidIndex][i] = p;
			}

			FORCE_INLINE const Real getKappa(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_kappa[fluidIndex][i];
			}

			FORCE_INLINE Real& getKappa(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_kappa[fluidIndex][i];
			}

			FORCE_INLINE void setKappa(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_kappa[fluidIndex][i] = p;
			}

			FORCE_INLINE const Real getKappaV(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_kappaV[fluidIndex][i];
			}

			FORCE_INLINE Real& getKappaV(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_kappaV[fluidIndex][i];
			}

			FORCE_INLINE void setKappaV(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_kappaV[fluidIndex][i] = p;
			}

			FORCE_INLINE const Real getDensityAdv(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_density_adv[fluidIndex][i];
			}

			FORCE_INLINE Real& getDensityAdv(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_density_adv[fluidIndex][i];
			}

			FORCE_INLINE void setDensityAdv(const unsigned int fluidIndex, const unsigned int i, const Real d)
			{
				m_density_adv[fluidIndex][i] = d;
			}

			FORCE_INLINE const Vector3r& getOldVelocity(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_oldVelocity[fluidIndex][i];
			}

			FORCE_INLINE void setOldVelocity(const unsigned int fluidIndex, const unsigned int i, const Vector3r& velocity)
			{
				m_oldVelocity[fluidIndex][i] = velocity;
			}

			FORCE_INLINE Vector3r& getAvgAcceleration(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_averageAcceleration[fluidIndex][i];
			}

			FORCE_INLINE Vector3r& getAvgAccelerationBorder(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_averageAccelerationBorder[fluidIndex][i];
			}
	};
}

#endif
