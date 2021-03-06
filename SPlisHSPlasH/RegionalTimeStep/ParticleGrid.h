#ifndef __ParticleGrid_h__
#define __ParticleGrid_h__

#include "SPlisHSPlasH/Common.h"

#include <array>
#include <vector>

#define REGION_LEVELS_COUNT 3
// Equivalent to the 'N' variable in Koike et al. Except, defining regional levels is done using the multiplier below.
#define LEVEL_TIMESTEP_MULTIPLIER 2
#define REGION_SPEED_MULTIPLIER 1.45f

namespace SPH
{
    struct CellParticlePairIndices
    {
        /* Indices for m_cellParticlePairs. These mark the interval for the cell's particles within said array. */
        unsigned int pairIndexBegin;
        unsigned int pairIndexEnd;
    };

    class ParticleGrid
    {
        public:
            ParticleGrid() = default;
            ~ParticleGrid();

            void init();

            /** Places particles into cells, and calculates the max velocity for each cell.
             * @param levelMinSpeeds The minimum speed for each region level. Particle's speed must be faster than the
             * minimum to belong to a level. The first level is the fastest one.
            */
            void determineRegions();

            // Writes border particle indices to m_particleIndices
            void enableBorderParticleIndices(unsigned int modelIdx, unsigned int level);

            void defineLevelParticleIndices();

            // Defines particle indices for a level and its sublevels
            void defineLevelParticleIndices(unsigned int highestLevel);

            FORCE_INLINE const std::vector<unsigned int> *getLevelParticleIndices() const           { return m_particleIndices.data(); }
            FORCE_INLINE const unsigned int *getLevelParticleCounts(unsigned int level) const       { return m_levelParticleCounts[level].data(); }
            FORCE_INLINE const unsigned int *getLevelUnionParticleCounts(unsigned int level) const  { return m_levelUnionsParticleCounts[level].data(); }
            FORCE_INLINE unsigned int getLevelBorderParticleCounts(unsigned int level, unsigned int modelIdx) const  { return m_levelBorderParticleCounts[level][modelIdx]; }
            FORCE_INLINE const Vector3f& getGridSize() const { return m_gridSize; }
            FORCE_INLINE unsigned int getParticleLevel(unsigned int modelIdx, unsigned int particleIdx) const { return m_particleLevels[modelIdx][particleIdx]; }

            FORCE_INLINE bool isBorder(unsigned int fluidModelIndex, unsigned int i) const { return m_particleBorderLevels[fluidModelIndex][i] != UINT32_MAX; }
            FORCE_INLINE const unsigned int* getParticleBorderLevels(unsigned int modelIdx) const { return m_particleBorderLevels[modelIdx].data(); }

        private:
            void initGridSizeAndResolution();

            void resetCells();

            // Puts particles into cells and finds the max velocity in each cell
            void findCellParticlePairs();

            // Iterates through all grid cells to find the largest velocity. Results are also stored in m_maxSpeedSquared.
            void findMaxSpeedSquared();

            // Should only be called before the cell-particle pairs array is sorted
            void defineCellLevels();
            void defineParticleLevels();

            void identifyRegionBorders();

            // Writes border particle indices to m_borderParticleIndices
            void writeBorderIndices();

        private:
            Vector3i m_resolution;

            struct CellParticlePair
            {
                unsigned int cellIndex;
                unsigned int particleIndex;

                bool operator<(const CellParticlePair& other)
                {
                    return cellIndex < other.cellIndex;
                }
            };

            // Each particle's cell index
            std::vector<std::vector<CellParticlePair>> m_cellParticlePairs;
            std::vector<std::vector<unsigned int>> m_particleLevels;

            // Particle indices, sorted by regional level
            std::vector<std::vector<unsigned int>> m_particleIndices;

            /*  Indices to particles that belong to borders. Sorted by the lowest neighboring level.
                Use m_levelBorderParticleCounts to find each level's bordering particles. */
            std::vector<std::vector<unsigned int>> m_borderParticleIndices;

            /*  The amount of particles per level. The latter contains the particle count of each level and their
                sublevels. */
            std::array<std::vector<unsigned int>, REGION_LEVELS_COUNT> m_levelParticleCounts;
            std::array<std::vector<unsigned int>, REGION_LEVELS_COUNT> m_levelUnionsParticleCounts;

            std::vector<std::vector<Real>> m_cellSpeedsSquared;
            std::vector<std::vector<unsigned int>> m_cellRegionLevels;
            std::vector<std::vector<CellParticlePairIndices>> m_cellPairIndices;

            /*  One element per cell. If a cell is in a region border, its regional level will be stored. If the cell
                is not in a border, UINT32_MAX will be stored instead. */
            std::vector<std::vector<unsigned int>> m_cellBorderLevels;

            // One element per particle. 1 signifies that a particle is in a regional border
            std::vector<std::vector<unsigned int>> m_particleBorderLevels;

            std::array<std::vector<unsigned int>, REGION_LEVELS_COUNT - 1> m_levelBorderParticleCounts;

            // Assumes the scene uses unitbox.obj. Each component is a distance to a wall.
            Vector3r m_gridSize = { 1.0f, 1.0f, 1.0f };

            // The maximum velocity out of every fluid particle
            Real m_maxSpeedSquared = -1.0f;
    };
}

#endif
