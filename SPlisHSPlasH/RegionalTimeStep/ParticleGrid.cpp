#include "ParticleGrid.h"

#include "Simulator/SceneConfiguration.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Simulator/SimulatorBase.h"

using namespace SPH;

void ParticleGrid::init()
{
    Utilities::SceneLoader::Scene &scene = SceneConfiguration::getCurrent()->getScene();
    if (!scene.useRegionalTimeStepping)
        return;

    initGridSizeAndResolution();

    Simulation *simulation = Simulation::getCurrent();

    const unsigned int nModels = simulation->numberOfFluidModels();
    m_cellParticlePairs.resize(nModels);
    for (unsigned int modelIndex = 0; modelIndex < nModels; modelIndex++)
    {
        const FluidModel *fluidModel = simulation->getFluidModel(modelIndex);
		m_cellParticlePairs[modelIndex].resize(fluidModel->numActiveParticles());
    }
}

void ParticleGrid::determineRegions()
{
    resetCells();
    findCellParticlePairs();
    findMaxVelocity();

    const Real smallestTimeStep = TimeManager::getCurrent()->getTimeStepSize();

    // Decide each cell's level
    for (std::vector<GridCell> &fluidModelCells : m_cells)
    {
        #pragma omp parallel default(shared)
        {
            // Use each block's max velocity to determine its level
            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < fluidModelCells.size(); cellIdx++)
            {
                GridCell &cell = fluidModelCells[cellIdx];
                cell.regionLevel = unsigned int(m_maxVelocity / cell.maxSpeed) / LEVEL_TIMESTEP_MULTIPLIER;
            }
        }
    }
}

void ParticleGrid::simulateLevels()
{
    m_levels[0].substepsRemaining = 1;

    for (unsigned int levelIdx = 1; levelIdx < REGION_LEVELS_COUNT; levelIdx++)
    {
        m_levels[levelIdx].substepsRemaining = m_levels[levelIdx - 1].substepsRemaining * LEVEL_TIMESTEP_MULTIPLIER;
    }

    // simulateLevel(REGION_LEVELS_COUNT - 1);
    // std::vector<unsigned int> interp
}

void ParticleGrid::initGridSizeAndResolution()
{
    Utilities::SceneLoader::Scene &scene = SceneConfiguration::getCurrent()->getScene();

    // Find wall distances
    m_gridSize = { 1.0f, 1.0f, 1.0f };
    for (Utilities::SceneLoader::BoundaryData* boundaryModel : scene.boundaryModels)
    {
        if (boundaryModel->isWall)
        {
            m_gridSize = boundaryModel->scale;
            break;
        }
    }

    Simulation *simulation = Simulation::getCurrent();
    const Real cellSize = simulation->getSupportRadius();

    m_resolution = (m_gridSize * (1.0f / cellSize)).cast<int>();

    const int cellCount = m_resolution.x() * m_resolution.y() * m_resolution.z();
    m_cells.resize(simulation->numberOfFluidModels());
    for (std::vector<GridCell> &fluidModelGrid : m_cells)
    {
        fluidModelGrid.resize(cellCount);
    }
}

void ParticleGrid::resetCells()
{
    // Reset each cell's max speed
    for (std::vector<GridCell> &fluidModelCells : m_cells)
    {
        #pragma omp parallel default(shared)
        {
            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < fluidModelCells.size(); cellIdx++)
            {
                GridCell &cell = fluidModelCells[cellIdx];
                cell.maxSpeed = 0.0f;
            }
        }
    }
}

void ParticleGrid::findCellParticlePairs()
{
    Simulation *simulation = Simulation::getCurrent();
    const int nModels = simulation->numberOfFluidModels();

    /*  Place particles into cells and define the max velocity for each cell.
        Placing particles into cells is a 3-step algorithm:
        1. Each particle is assigned a cell index. These cell indices are stored alongside particle indices as pairs
        in an array. Each cell's max speed is also calculated in this step.
        2. Sort the cell index array.
        3. Find each cell's interval of particle indices in the sorted array. */

    // Step 1
    for (int modelIndex = 0; modelIndex < nModels; modelIndex++)
	{
		#pragma omp parallel default(shared)
		{
		    const FluidModel *fluidModel = simulation->getFluidModel(modelIndex);
		    const int numParticles = fluidModel->numActiveParticles();

            // Reciprocal of twice the grid's dimensions
            const Vector3r gridDoubleSizeReciprocal = (m_gridSize * 2.0f).cwiseInverse();
            const Real timeStepSize = TimeManager::getCurrent()->getTimeStepSize();

            std::vector<GridCell> &fluidModelCells = m_cells[modelIndex];
            std::vector<CellParticlePair>& cellParticleIndices = m_cellParticlePairs[modelIndex];

            const Vector3r gridResReal = m_resolution.cast<Real>();

			#pragma omp for schedule(static)
			for (int i = 0; i < numParticles; i++)
			{
				if (fluidModel->getParticleState(i) == ParticleState::Active)
				{
                    // Place the particle in a cell
                    const Vector3r &particlePosition = fluidModel->getPosition(i);

                    // The grid size is added to the particle's position to avoid negaive numbers
                    const Vector3i indices = (particlePosition + m_gridSize).cwiseProduct(gridDoubleSizeReciprocal).cwiseProduct(gridResReal).cast<int>();

                    const unsigned int cellIndex =
                        indices.x() +
                        indices.y() * m_resolution.x() +
                        indices.z() * m_resolution.x() * m_resolution.y();

                    cellParticleIndices[i] = { cellIndex, (unsigned int)i };

                    // Update cell's max speed
                    GridCell &cell = fluidModelCells[cellIndex];
					const Vector3r &particleVelocity = fluidModel->getVelocity(i);
                    const Vector3r &particleAcceleration = fluidModel->getAcceleration(i);
                    const Real particleSpeed = (particleVelocity + particleAcceleration * timeStepSize).squaredNorm();

                    cell.maxSpeed = std::max(particleSpeed, cell.maxSpeed);
				}
			}
        }
    }

    // Step 2
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)
        for (int modelIndex = 0; modelIndex < nModels; modelIndex++)
	    {
            std::vector<CellParticlePair> &cellParticlePairs = m_cellParticlePairs[modelIndex];
            std::sort(cellParticlePairs.begin(), cellParticlePairs.end());
        }
    }

    // Step 3
    for (int modelIndex = 0; modelIndex < nModels; modelIndex++)
    {
        std::vector<GridCell> &gridCells = m_cells[modelIndex];
        const std::vector<CellParticlePair> &cellParticlePairs = m_cellParticlePairs[modelIndex];

        unsigned int cellParticlePairIndex = 0;

        for (unsigned int cellIndex = 0; cellIndex < gridCells.size(); cellIndex++)
        {
            GridCell &cell = gridCells[cellIndex];
            cell.cellIndexBegin = cellParticlePairIndex;

            const unsigned int beginPairIndex = cellParticlePairIndex;
            while (cellParticlePairIndex < cellParticlePairs.size() && cellParticlePairs[cellParticlePairIndex].cellIndex == cellIndex)
                cellParticlePairIndex += 1;

            cell.cellIndexEnd = cellParticlePairIndex;
        }
    }
}

Real ParticleGrid::findMaxVelocity()
{
    for (const std::vector<GridCell> &fluidModelCells : m_cells)
    {
        for (const GridCell &cell : fluidModelCells)
        {
            m_maxVelocity = std::max(m_maxVelocity, cell.maxSpeed);
        }
    }

    return m_maxVelocity;
}

void ParticleGrid::simulateLevel(unsigned int level)
{
    // std::vector<unsigned int> interp
}
