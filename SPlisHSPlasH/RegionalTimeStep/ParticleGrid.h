#ifndef __ParticleGrid_h__
#define __ParticleGrid_h__

#include "SPlisHSPlasH/Common.h"

#include <array>
#include <vector>

#define REGION_LEVELS_COUNT 3
// Equivalent to the 'N' variable in Koike et al.
#define LEVEL_TIMESTEP_MULTIPLIER 2

namespace SPH
{
    struct Level
    {
        Real currentTime;
        Real nextTimestep;
        uint32_t substepsRemaining; // Substeps remaining in the current timestep
        std::vector<unsigned int> regions;
    };

    struct Region
    {
        std::vector<unsigned int> gridCells;
        std::vector<unsigned int> borderingRegions;
    };

    struct GridCell
    {
        Real maxSpeed;
        unsigned int regionLevel;

        /* Indices for m_cellParticlePairs. These mark the interval for the cell's particles within said array. */
        unsigned int pairIndexBegin;
        unsigned int pairIndexEnd;
    };

    class ParticleGrid
    {
        public:
            ParticleGrid() = default;
            ~ParticleGrid() = default;

            void init();

            /** Places particles into cells, and calculates the max velocity for each cell.
             * @param levelMinSpeeds The minimum speed for each region level. Particle's speed must be faster than the
             * minimum to belong to a level. The first level is the fastest one.
            */
            void determineRegions();

            void simulateLevels();

            void toggleRegionColorsRendering(bool enabled);

            FORCE_INLINE Real getMaxVelocity() const        { return m_maxVelocity; }
            FORCE_INLINE bool regionColorsEnabled()    { return m_regionColorsEnabled; }

        private:
            void initGridSizeAndResolution();

            void resetCells();

            // Puts particles into cells and finds the max velocity in each cell
            void findCellParticlePairs();

            // Iterates through all grid cells to find the largest velocity. Results are also stored in m_maxVelocity.
            Real findMaxVelocity();

            // Recursively simulates each level in the simulation
            void simulateLevel(unsigned int level);

            // Should only be called before the cell-particle pairs array is sorted
            void defineCellLevels();
            void defineParticleLevels();

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
            // Each particle's level. This is only used for rendering particle region colors, and can be toggled.
            std::vector<std::vector<unsigned int>> m_particleLevels;

            std::vector<std::vector<GridCell>> m_cells;

            // Assumes the scene uses unitbox.obj. Each component is a distance to a wall.
            Vector3r m_gridSize = { 1.0f, 1.0f, 1.0f };

            std::array<Level, REGION_LEVELS_COUNT> m_levels;
            std::vector<Region> m_regions;

            // The maximum velocity out of every fluid particle
            Real m_maxVelocity = -1.0f;

            bool m_regionColorsEnabled = false;
    };
}

#endif
