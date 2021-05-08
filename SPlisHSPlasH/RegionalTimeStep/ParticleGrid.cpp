#include "ParticleGrid.h"

#include "Simulator/SceneConfiguration.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Simulator/SimulatorBase.h"

using namespace SPH;

ParticleGrid::~ParticleGrid()
{
    Simulation *simulation = Simulation::getCurrent();
    const unsigned int nModels = simulation->numberOfFluidModels();

    for (unsigned int modelIdx = 0; modelIdx < nModels; modelIdx++)
    {
        FluidModel *fluidModel = simulation->getFluidModel(modelIdx);
        fluidModel->removeFieldByName("particle levels");

        if (!m_particleBorderLevels.empty())
            fluidModel->removeFieldByName("particle border levels");
    }
}

void ParticleGrid::init()
{
    initGridSizeAndResolution();

    Simulation *simulation = Simulation::getCurrent();

    const unsigned int nModels = simulation->numberOfFluidModels();
    m_cellParticlePairs.resize(nModels);
    m_cellBorderLevels.resize(nModels);
    m_particleIndices.resize(nModels);
    m_borderParticleIndices.resize(nModels);
    m_particleLevels.resize(nModels);
    m_particleBorderLevels.resize(nModels);

    for (unsigned int modelIndex = 0; modelIndex < nModels; modelIndex++)
    {
        FluidModel *fluidModel = simulation->getFluidModel(modelIndex);
        const unsigned int nParticles = fluidModel->numActiveParticles();
		m_cellParticlePairs[modelIndex].resize(nParticles);
        m_cellBorderLevels[modelIndex].resize(m_cells[modelIndex].size());
        m_particleIndices[modelIndex].resize(nParticles);
        m_borderParticleIndices[modelIndex].reserve(nParticles);

        m_particleLevels[modelIndex].resize(nParticles);
        m_particleBorderLevels[modelIndex].resize(nParticles);
        fluidModel->addField({ "particle levels", FieldType::UInt, [this, modelIndex](unsigned int particleIdx) { return &m_particleLevels[modelIndex][particleIdx]; } });
        fluidModel->addField({ "particle border levels", FieldType::UInt, [this, modelIndex](unsigned int particleIdx) { return &m_particleBorderLevels[modelIndex][particleIdx]; } });
    }

    for (unsigned int level = 0; level < REGION_LEVELS_COUNT; level++)
    {
        m_levelParticleCounts[level].resize(nModels);
        m_levelUnionsParticleCounts[level].resize(nModels);
    }

    for (unsigned int level = 0; level < REGION_LEVELS_COUNT - 1; level++)
        m_levelBorderParticleCounts[level].resize(nModels);
}

void ParticleGrid::determineRegions()
{
    resetCells();
    findCellParticlePairs();
    findMaxSpeedSquared();
    defineCellLevels();
    defineParticleLevels();

    /* These two functions could be executed in parallel */
    identifyRegionBorders();
    defineLevelParticleIndices();

    writeBorderIndices();
}

void ParticleGrid::calculateLevelBorder(unsigned int level)
{
    Simulation* sim = Simulation::getCurrent();
    const unsigned int nModels = sim->numberOfFluidModels();

    for (unsigned int modelIdx = 0; modelIdx < nModels; modelIdx++)
    {
        const std::vector<GridCell>& cells = m_cells[modelIdx];
        const std::vector<CellParticlePair>& cellParticlePairs = m_cellParticlePairs[modelIdx];

        // Start writing border particle indices where the higher level indices are currently stored
        const unsigned int indexStartPos = m_levelUnionsParticleCounts[level][modelIdx];
        unsigned int borderParticleCount = 0;

        if (indexStartPos < sim->getFluidModel(modelIdx)->getNumActiveParticles0())
        {
            unsigned int* borderParticleIndices = &m_particleIndices[modelIdx][indexStartPos];

            const std::vector<unsigned int>& borderLevels = m_cellBorderLevels[modelIdx];

            for (unsigned int cellIdx = 0; cellIdx < cells.size(); cellIdx++)
            {
                const unsigned int borderLevel = borderLevels[cellIdx];
                if (borderLevel == level)
                {
                    const GridCell& cell = cells[cellIdx];
                    unsigned int pairIdx = cell.pairIndexBegin;
                    const unsigned int pairIdxEnd = cell.pairIndexEnd;

                    while (pairIdx < pairIdxEnd)
                    {
                        const unsigned int particleIdx = cellParticlePairs[pairIdx++].particleIndex;
                        borderParticleIndices[borderParticleCount++] = particleIdx;
                    }
                }
            }
        }

        // Update particle counts
        m_levelBorderParticleCounts[level][modelIdx] = borderParticleCount;
    }
}

void ParticleGrid::enableBorderParticleIndices(unsigned int modelIdx, unsigned int level)
{
    // Start writing border particle indices where the higher level indices are currently stored
    const unsigned int indexStartPos = m_levelUnionsParticleCounts[level][modelIdx];
    const unsigned int levelBorderParticleCount = m_levelBorderParticleCounts[level][modelIdx];

    if (levelBorderParticleCount != 0 && indexStartPos < m_particleIndices[modelIdx].size())
    {
        // Calculate the index from where to start copying border particle indices
        unsigned int borderParticlesStartIdx = 0;
        for (unsigned int lowerLevel = 0; lowerLevel < level; lowerLevel++)
        {
            borderParticlesStartIdx += m_levelBorderParticleCounts[lowerLevel][modelIdx];
        }

        unsigned int* particleIndices = &m_particleIndices[modelIdx][indexStartPos];
        const unsigned int* borderParticleIndices = &m_borderParticleIndices[modelIdx][borderParticlesStartIdx];

        std::copy_n(borderParticleIndices, levelBorderParticleCount, particleIndices);
    }
}

void ParticleGrid::defineLevelParticleIndices()
{
    Simulation *simulation = Simulation::getCurrent();
    const unsigned int numFluidModels = simulation->numberOfFluidModels();

    for (unsigned int modelIdx = 0; modelIdx < numFluidModels; modelIdx++)
    {
        // #pragma omp parallel num_threads(2) default(shared)
        {
            const std::vector<unsigned int> &particleLevels = m_particleLevels[modelIdx];
            std::vector<unsigned int>& levelParticleIndices = m_particleIndices[modelIdx];

            unsigned int levelIndicesPos = 0;

            // const int threadID = omp_get_thread_num();
            // if (threadID == 0)
            for (unsigned int level = 0; level < REGION_LEVELS_COUNT; level++)
            {
                const unsigned int indicesPosBegin = levelIndicesPos;

                const unsigned int particleCount = particleLevels.size();
                for (unsigned int particleIdx = 0; particleIdx < particleCount; particleIdx++)
                {
                    if (particleLevels[particleIdx] == level)
                    {
                        levelParticleIndices[levelIndicesPos++] = particleIdx;
                    }
                }

                m_levelParticleCounts[level][modelIdx] = levelIndicesPos - indicesPosBegin;
                m_levelUnionsParticleCounts[level][modelIdx] = m_levelParticleCounts[level][modelIdx] +
                    (level > 0 ? m_levelUnionsParticleCounts[level - 1][modelIdx] : 0);
            }
        }
    }
}

void ParticleGrid::defineLevelParticleIndices(unsigned int highestLevel)
{
    Simulation *simulation = Simulation::getCurrent();
    const unsigned int numFluidModels = simulation->numberOfFluidModels();

    for (unsigned int level = 1; level <= highestLevel; level++)
    {
        for (unsigned int modelIdx = 0; modelIdx < numFluidModels; modelIdx++)
        {
            // #pragma omp parallel num_threads(2) default(shared)
            const std::vector<unsigned int> &particleLevels = m_particleLevels[modelIdx];
            const unsigned int startingIndex = level > 0 ? m_levelUnionsParticleCounts[level - 1][modelIdx] : 0;
            if (startingIndex >= m_particleIndices[modelIdx].size())
            {
                continue;
            }

            unsigned int* levelParticleIndices = &m_particleIndices[modelIdx][startingIndex];

            unsigned int levelIndicesPos = 0;

            // const int threadID = omp_get_thread_num();
            // if (threadID == 0)
            std::vector<unsigned int> &levelUnionsParticleCounts = m_levelUnionsParticleCounts[level];

            const unsigned int indicesPosBegin = levelIndicesPos;

            const unsigned int particleCount = particleLevels.size();
            for (unsigned int particleIdx = 0; particleIdx < particleCount; particleIdx++)
            {
                if (particleLevels[particleIdx] == level)
                {
                    levelParticleIndices[levelIndicesPos++] = particleIdx;
                }
            }

            m_levelParticleCounts[level][modelIdx] = levelIndicesPos - indicesPosBegin;
            m_levelUnionsParticleCounts[level][modelIdx] = m_levelParticleCounts[level][modelIdx] +
                (level > 0 ? m_levelUnionsParticleCounts[level - 1][modelIdx] : 0);
        }
    }
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
                cell.maxSpeedSquared = 0.0f;
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
                // TODO: Move this, data race!
                GridCell &cell = fluidModelCells[cellIndex];
                const Vector3r &particleVelocity = fluidModel->getVelocity(i);
                const Vector3r &particleAcceleration = fluidModel->getAcceleration(i);
                const Real particleSpeedSquared = (particleVelocity + particleAcceleration * timeStepSize).squaredNorm();

                cell.maxSpeedSquared = std::max(particleSpeedSquared, cell.maxSpeedSquared);
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
            cell.pairIndexBegin = cellParticlePairIndex;

            while (cellParticlePairIndex < cellParticlePairs.size() && cellParticlePairs[cellParticlePairIndex].cellIndex == cellIndex)
                cellParticlePairIndex += 1;

            cell.pairIndexEnd = cellParticlePairIndex;
        }
    }
}

Real ParticleGrid::findMaxSpeedSquared()
{
    m_maxSpeedSquared = -1.0f;

    for (const std::vector<GridCell> &fluidModelCells : m_cells)
    {
        for (const GridCell &cell : fluidModelCells)
        {
            m_maxSpeedSquared = std::max(m_maxSpeedSquared, cell.maxSpeedSquared);
        }
    }

    // printf("Max Speed: %f\n", std::sqrtf(m_maxSpeedSquared));
    return m_maxSpeedSquared;
}

void ParticleGrid::defineCellLevels()
{
    for (std::vector<GridCell> &fluidModelCells : m_cells)
    {
        #pragma omp parallel default(shared)
        {
            const Real logN = std::logf(REGION_SPEED_MULTIPLIER);
            const Real timeStepSize = TimeManager::getCurrent()->getTimeStepSize();
            const Real maxTimeStepSize = Simulation::getCurrent()->getMaxTimeStepSize();

            // Limit what regional levels are enabled this frame, to stay within the max time step size
            const unsigned int maxRegionalLevel = std::min<unsigned int>(REGION_LEVELS_COUNT - 1, unsigned int(std::logf(maxTimeStepSize / timeStepSize) / logN));

            const Real maxSpeed = std::sqrtf(m_maxSpeedSquared);

            // Use each block's max velocity to determine its level
            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < fluidModelCells.size(); cellIdx++)
            {
                GridCell &cell = fluidModelCells[cellIdx];
                // cell.regionLevel = cellIdx % 3;
                cell.regionLevel = unsigned int(std::logf(maxSpeed / std::sqrtf(cell.maxSpeedSquared)) / logN);
                cell.regionLevel = std::min(cell.regionLevel, maxRegionalLevel);
            }
        }
    }
}

void ParticleGrid::defineParticleLevels()
{
    Simulation *simulation = Simulation::getCurrent();
    const unsigned int numFluidModels = simulation->numberOfFluidModels();

    for (unsigned int modelIdx = 0; modelIdx < numFluidModels; modelIdx++)
    {
        #pragma omp parallel default(shared)
        {
            std::vector<unsigned int> &particleLevels = m_particleLevels[modelIdx];
            const std::vector<GridCell> &gridCells = m_cells[modelIdx];
            const std::vector<CellParticlePair> &cellParticlePairs = m_cellParticlePairs[modelIdx];

            const int cellCount = (int)gridCells.size();

            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < cellCount; cellIdx++)
            {
                const GridCell &cell = gridCells[cellIdx];
                unsigned int pairIdx = cell.pairIndexBegin;
                const unsigned int pairIdxEnd = cell.pairIndexEnd;

                while (pairIdx < pairIdxEnd)
                {
                    const unsigned int particleIdx = cellParticlePairs[pairIdx++].particleIndex;
                    particleLevels[particleIdx] = cell.regionLevel;
                }
            }
        }
    }
}

void ParticleGrid::identifyRegionBorders()
{
    Simulation *simulation = Simulation::getCurrent();
    const unsigned int numFluidModels = simulation->numberOfFluidModels();

    for (unsigned int modelIdx = 0; modelIdx < numFluidModels; modelIdx++)
    {
        std::vector<unsigned int> &cellBorderLevels = m_cellBorderLevels[modelIdx];
        std::fill_n(cellBorderLevels.data(), cellBorderLevels.size(), UINT32_MAX);

        std::vector<unsigned int>& isBorder = m_particleBorderLevels[modelIdx];
        std::fill_n(isBorder.data(), isBorder.size(), UINT32_MAX);
    }

    for (unsigned int modelIdx = 0; modelIdx < numFluidModels; modelIdx++)
    {
        #pragma omp parallel default(shared)
        {
            std::vector<unsigned int> &cellBorderLevels = m_cellBorderLevels[modelIdx];
            std::vector<GridCell> &gridCells = m_cells[modelIdx];
            const int cellCount = (int)gridCells.size();
            const int zStep = m_resolution.x() * m_resolution.y();

            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < cellCount; cellIdx++)
            {
                const GridCell &cell = gridCells[cellIdx];

                // Skip empty cells
                if (cell.pairIndexBegin == cell.pairIndexEnd)
                    continue;

                const int cellPosZ = cellIdx / zStep;
                const int cellPosY = (cellIdx - cellPosZ * zStep) / m_resolution.x();
                const int cellPosX = cellIdx % m_resolution.x();

                const int offsetStartX = std::max(-1, -cellPosX);
                const int offsetStartY = std::max(-1, -cellPosY);
                const int offsetStartZ = std::max(-1, -cellPosZ);

                const int offsetEndX = std::min(2, m_resolution.x() - cellPosX);
                const int offsetEndY = std::min(2, m_resolution.y() - cellPosY);
                const int offsetEndZ = std::min(2, m_resolution.z() - cellPosZ);

                const unsigned int cellLevel = cell.regionLevel;

                // Check all neighbors in a 3x3x3 grid around the cell
                for (int offsetZ = offsetStartZ; offsetZ < offsetEndZ; offsetZ++)
                {
                    for (int offsetY = offsetStartY; offsetY < offsetEndY; offsetY++)
                    {
                        for (int offsetX = offsetStartX; offsetX < offsetEndX; offsetX++)
                        {
                            const unsigned int neighborIdx = cellIdx + offsetX + offsetY * m_resolution.x() + offsetZ * zStep;
                            const GridCell &neighborCell = gridCells[neighborIdx];

                            /*  If the neighboring cell is not empty (of particles) and its level is different, then
                                this cell is part of a border. */
                            if (neighborCell.pairIndexBegin != neighborCell.pairIndexEnd && neighborCell.regionLevel < cellLevel)
                            {
                                cellBorderLevels[cellIdx] = std::min(cellBorderLevels[cellIdx], neighborCell.regionLevel);
                            }
                        }
                    }
                }

                const unsigned int cellBorderLevel = cellBorderLevels[cellIdx];
                if (cellBorderLevel != UINT32_MAX)
                {
                    unsigned int pairIdx = cell.pairIndexBegin;
                    const unsigned int pairIdxEnd = cell.pairIndexEnd;

                    std::vector<unsigned int>& particleBorderLevels = m_particleBorderLevels[modelIdx];
                    const std::vector<CellParticlePair>& cellParticlePairs = m_cellParticlePairs[modelIdx];

                    while (pairIdx < pairIdxEnd)
                    {
                        const unsigned int particleIdx = cellParticlePairs[pairIdx++].particleIndex;
                        particleBorderLevels[particleIdx] = cellBorderLevel;
                    }
                }
            }
        }
    }
}

void ParticleGrid::writeBorderIndices()
{
    Simulation *simulation = Simulation::getCurrent();
    const unsigned int numFluidModels = simulation->numberOfFluidModels();

    for (unsigned int modelIdx = 0; modelIdx < numFluidModels; modelIdx++)
    {
        std::vector<unsigned int>& borderParticleIndices = m_borderParticleIndices[modelIdx];
        borderParticleIndices.clear();

        const std::vector<unsigned int> &borderLevels = m_cellBorderLevels[modelIdx];
        const std::vector<GridCell> &gridCells = m_cells[modelIdx];
        const std::vector<CellParticlePair>& cellParticlePairs = m_cellParticlePairs[modelIdx];

        const unsigned int cellCount = borderLevels.size();

        for (unsigned int level = 0; level < REGION_LEVELS_COUNT - 1; level++)
        {
            const unsigned int indicesInitialSize = borderParticleIndices.size();

            for (unsigned int cellIdx = 0; cellIdx < cellCount; cellIdx++)
            {
                if (borderLevels[cellIdx] == level)
                {
                    const GridCell &cell = gridCells[cellIdx];

                    unsigned int pairIdx = cell.pairIndexBegin;
                    const unsigned int pairIdxEnd = cell.pairIndexEnd;

                    while (pairIdx < pairIdxEnd)
                    {
                        const unsigned int particleIdx = cellParticlePairs[pairIdx++].particleIndex;
                        borderParticleIndices.push_back(particleIdx);
                    }
                }
            }

            m_levelBorderParticleCounts[level][modelIdx] = borderParticleIndices.size() - indicesInitialSize;
        }
    }
}
