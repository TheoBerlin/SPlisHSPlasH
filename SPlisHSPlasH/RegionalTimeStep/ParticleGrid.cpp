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
    const int cellCount = m_resolution.x() * m_resolution.y() * m_resolution.z();

    Simulation *simulation = Simulation::getCurrent();

    const unsigned int nModels = simulation->numberOfFluidModels();
    m_cellParticlePairs.resize(nModels);
    m_cellBorderLevels.resize(nModels);
    m_cellSecondBorderLevels.resize(nModels);
    m_particleIndices.resize(nModels);
    m_borderParticleIndices.resize(nModels);
    m_secondBorderParticleIndices.resize(nModels);
    m_particleLevels.resize(nModels);
    m_particleBorderLevels.resize(nModels);
    m_cellSpeedsSquared.resize(nModels);
    m_cellRegionLevels.resize(nModels);
    m_cellPairIndices.resize(nModels);

    for (unsigned int modelIdx = 0; modelIdx < nModels; modelIdx++)
    {
        FluidModel *fluidModel = simulation->getFluidModel(modelIdx);
        const unsigned int nParticles = fluidModel->numActiveParticles();
		m_cellParticlePairs[modelIdx].resize(nParticles);
        m_particleIndices[modelIdx].resize(nParticles);
        m_borderParticleIndices[modelIdx].reserve(nParticles);
        m_secondBorderParticleIndices[modelIdx].reserve(nParticles);

        m_cellSpeedsSquared[modelIdx].resize(cellCount);
        m_cellRegionLevels[modelIdx].resize(cellCount);
        m_cellPairIndices[modelIdx].resize(cellCount);
        m_cellBorderLevels[modelIdx].resize(cellCount);
        m_cellSecondBorderLevels[modelIdx].resize(cellCount);

        m_particleLevels[modelIdx].resize(nParticles);
        m_particleBorderLevels[modelIdx].resize(nParticles);
        fluidModel->addField({ "particle levels", FieldType::UInt, [this, modelIdx](unsigned int particleIdx) { return &m_particleLevels[modelIdx][particleIdx]; } });
        fluidModel->addField({ "particle border levels", FieldType::UInt, [this, modelIdx](unsigned int particleIdx) { return &m_particleBorderLevels[modelIdx][particleIdx]; } });
    }

    for (unsigned int level = 0; level < REGION_LEVELS_COUNT; level++)
    {
        m_levelParticleCounts[level].resize(nModels);
        m_levelUnionsParticleCounts[level].resize(nModels);
    }

    for (unsigned int level = 0; level < REGION_LEVELS_COUNT - 1; level++)
    {
        m_levelBorderParticleCounts[level].resize(nModels);
        m_nSecondBorderUnionParticles[level].resize(nModels);
    }
}

void ParticleGrid::determineRegions()
{
    resetCells();
    findCellParticlePairs();
    findMaxSpeedSquared();
    defineCellLevels();
    defineParticleLevels();

    #pragma omp parallel sections
    {
        {
            identifyRegionBorders();
            writeBorderIndices();
        }

        #pragma omp section
        {
            defineLevelParticleIndices();
        }
    }
}

void ParticleGrid::calculateLevelBorder(unsigned int level)
{
    Simulation* sim = Simulation::getCurrent();
    const unsigned int nModels = sim->numberOfFluidModels();

    for (unsigned int modelIdx = 0; modelIdx < nModels; modelIdx++)
    {
        const std::vector<CellParticlePairIndices>& pairIndices = m_cellPairIndices[modelIdx];
        const std::vector<CellParticlePair>& cellParticlePairs = m_cellParticlePairs[modelIdx];

        // Start writing border particle indices where the higher level indices are currently stored
        const unsigned int indexStartPos = m_levelUnionsParticleCounts[level][modelIdx];
        unsigned int borderParticleCount = 0;

        if (indexStartPos < sim->getFluidModel(modelIdx)->getNumActiveParticles0())
        {
            unsigned int* borderParticleIndices = &m_particleIndices[modelIdx][indexStartPos];

            const std::vector<unsigned int>& borderLevels = m_cellBorderLevels[modelIdx];

            for (unsigned int cellIdx = 0; cellIdx < pairIndices.size(); cellIdx++)
            {
                const unsigned int borderLevel = borderLevels[cellIdx];
                if (borderLevel == level)
                {
                    const CellParticlePairIndices& cell = pairIndices[cellIdx];
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

        // Write second border particle indices
        const unsigned int* secondBorderParticles = m_secondBorderParticleIndices[modelIdx].data();
        const unsigned int secondBorderParticleCount = m_nSecondBorderUnionParticles[level][modelIdx];
        std::copy_n(secondBorderParticles, secondBorderParticleCount, &particleIndices[levelBorderParticleCount]);
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
}

void ParticleGrid::resetCells()
{
    const unsigned int nModels = m_cellPairIndices.size();

    // Reset each cell's max speed
    for (unsigned int modelIdx = 0; modelIdx < nModels; modelIdx++)
    {
        const unsigned int nCells = m_cellPairIndices[modelIdx].size();

        std::vector<Real>& cellSpeedsSquared = m_cellSpeedsSquared[modelIdx];
        std::fill_n(cellSpeedsSquared.data(), nCells, 0.0f);
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
    for (int modelIdx = 0; modelIdx < nModels; modelIdx++)
    {
        std::vector<CellParticlePairIndices>& cellParticlePairIndices = m_cellPairIndices[modelIdx];
        const std::vector<CellParticlePair>& cellParticlePairs = m_cellParticlePairs[modelIdx];

        unsigned int cellParticlePairIndex = 0;

        for (unsigned int cellIndex = 0; cellIndex < cellParticlePairIndices.size(); cellIndex++)
        {
            CellParticlePairIndices& pairIndices = cellParticlePairIndices[cellIndex];
            pairIndices.pairIndexBegin = cellParticlePairIndex;

            while (cellParticlePairIndex < cellParticlePairs.size() && cellParticlePairs[cellParticlePairIndex].cellIndex == cellIndex)
                cellParticlePairIndex += 1;

            pairIndices.pairIndexEnd = cellParticlePairIndex;
        }
    }
}

void ParticleGrid::findMaxSpeedSquared()
{
    Simulation* sim = Simulation::getCurrent();
    const unsigned int nModels = sim->numberOfFluidModels();

    // Update each cell's max speed
    for (unsigned int modelIdx = 0; modelIdx < nModels; modelIdx++)
    {
        #pragma omp parallel default(shared)
        {
            const FluidModel* model = sim->getFluidModel(modelIdx);
            const Real timeStepSize = TimeManager::getCurrent()->getTimeStepSize();

            std::vector<Real>& cellSpeedsSquared = m_cellSpeedsSquared[modelIdx];
            const std::vector<CellParticlePairIndices>& cellPairIndices = m_cellPairIndices[modelIdx];
            const std::vector<CellParticlePair>& cellParticlePairs = m_cellParticlePairs[modelIdx];

            const unsigned int nCells = cellSpeedsSquared.size();

            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < nCells; cellIdx++)
            {
                const CellParticlePairIndices& pairIndices = cellPairIndices[cellIdx];
                unsigned int pairIdx = pairIndices.pairIndexBegin;
                const unsigned int pairIdxEnd = pairIndices.pairIndexEnd;

                Real& cellMaxSpeedSquared = cellSpeedsSquared[cellIdx];

                while (pairIdx < pairIdxEnd)
                {
                    const unsigned int i = cellParticlePairs[pairIdx++].particleIndex;
                    const Vector3r &particleVelocity = model->getVelocity(i);
                    const Vector3r &particleAcceleration = model->getAcceleration(i);
                    const Real particleSpeedSquared = (particleVelocity + particleAcceleration * timeStepSize).squaredNorm();

                    cellMaxSpeedSquared = std::max(cellMaxSpeedSquared, particleSpeedSquared);
                }
            }
        }
    }

    m_maxSpeedSquared = -1.0f;

    for (const std::vector<Real>& cellSpeedsSquared : m_cellSpeedsSquared)
    {
        for (const Real cellMaxSpeedSquared : cellSpeedsSquared)
        {
            m_maxSpeedSquared = std::max(m_maxSpeedSquared, cellMaxSpeedSquared);
        }
    }

    // printf("Max Speed: %f\n", std::sqrtf(m_maxSpeedSquared));
}

void ParticleGrid::defineCellLevels()
{
    const unsigned int nModels = m_cellRegionLevels.size();

    for (unsigned int modelIdx = 0; modelIdx < nModels; modelIdx++)
    {
        #pragma omp parallel default(shared)
        {
            std::vector<unsigned int>& cellRegionLevels = m_cellRegionLevels[modelIdx];
            const std::vector<Real>& cellMaxSpeedsSquared = m_cellSpeedsSquared[modelIdx];

            const Real logN = std::logf(REGION_SPEED_MULTIPLIER);
            const Real timeStepSize = TimeManager::getCurrent()->getTimeStepSize();
            const Real maxTimeStepSize = Simulation::getCurrent()->getMaxTimeStepSize();

            // Limit what regional levels are enabled this frame, to stay within the max time step size
            const unsigned int maxRegionalLevel = std::min<unsigned int>(REGION_LEVELS_COUNT - 1, unsigned int(std::logf(maxTimeStepSize / timeStepSize) / logN));

            const Real maxSpeed = std::sqrtf(m_maxSpeedSquared);

            // Use each block's max velocity to determine its level
            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < cellRegionLevels.size(); cellIdx++)
            {
                unsigned int& regionLevel = cellRegionLevels[cellIdx];
                // cell.regionLevel = cellIdx % 3;
                regionLevel = unsigned int(std::logf(maxSpeed / std::sqrtf(cellMaxSpeedsSquared[cellIdx])) / logN);
                regionLevel = std::min(regionLevel, maxRegionalLevel);
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
            const std::vector<CellParticlePairIndices> &cellParticlePairIndices = m_cellPairIndices[modelIdx];
            const std::vector<CellParticlePair> &cellParticlePairs = m_cellParticlePairs[modelIdx];
            const std::vector<unsigned int> &cellRegionLevels = m_cellRegionLevels[modelIdx];

            const int cellCount = (int)cellParticlePairIndices.size();

            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < cellCount; cellIdx++)
            {
                const unsigned int regionLevel = cellRegionLevels[cellIdx];

                const CellParticlePairIndices& pairIndices = cellParticlePairIndices[cellIdx];
                unsigned int pairIdx = pairIndices.pairIndexBegin;
                const unsigned int pairIdxEnd = pairIndices.pairIndexEnd;

                while (pairIdx < pairIdxEnd)
                {
                    const unsigned int particleIdx = cellParticlePairs[pairIdx++].particleIndex;
                    particleLevels[particleIdx] = regionLevel;
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
            std::vector<unsigned int>& cellBorderLevels = m_cellBorderLevels[modelIdx];
            std::vector<unsigned int>& cellRegionLevels = m_cellRegionLevels[modelIdx];
            std::vector<CellParticlePairIndices>& cellParticlePairIndices = m_cellPairIndices[modelIdx];

            const int cellCount = (int)cellBorderLevels.size();
            const int zStep = m_resolution.x() * m_resolution.y();

            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < cellCount; cellIdx++)
            {
                const CellParticlePairIndices& pairIndices = cellParticlePairIndices[cellIdx];
                const unsigned int cellLevel = cellRegionLevels[cellIdx];

                // Skip empty cells and level 0 cells
                if (pairIndices.pairIndexBegin == pairIndices.pairIndexEnd || cellLevel == 0)
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

                // Check all neighbors in a 3x3x3 grid around the cell
                for (int offsetZ = offsetStartZ; offsetZ < offsetEndZ; offsetZ++)
                {
                    for (int offsetY = offsetStartY; offsetY < offsetEndY; offsetY++)
                    {
                        for (int offsetX = offsetStartX; offsetX < offsetEndX; offsetX++)
                        {
                            const unsigned int neighborIdx = cellIdx + offsetX + offsetY * m_resolution.x() + offsetZ * zStep;
                            const CellParticlePairIndices& neighborCellPairIndices = cellParticlePairIndices[neighborIdx];
                            const unsigned int neighborRegionLevel = cellRegionLevels[neighborIdx];

                            /*  If the neighboring cell is not empty (of particles) and its level is different, then
                                this cell is part of a border. */
                            if (neighborCellPairIndices.pairIndexBegin != neighborCellPairIndices.pairIndexEnd && neighborRegionLevel < cellLevel)
                            {
                                cellBorderLevels[cellIdx] = std::min(cellBorderLevels[cellIdx], neighborRegionLevel);

                                if (neighborRegionLevel == 0)
                                    goto assignParticleLevels;
                            }
                        }
                    }
                }

            assignParticleLevels: ;
                const unsigned int cellBorderLevel = cellBorderLevels[cellIdx];
                if (cellBorderLevel != UINT32_MAX)
                {
                    unsigned int pairIdx = pairIndices.pairIndexBegin;
                    const unsigned int pairIdxEnd = pairIndices.pairIndexEnd;

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

void ParticleGrid::identifySecondBorderCells()
{
    Simulation *simulation = Simulation::getCurrent();
    const unsigned int numFluidModels = simulation->numberOfFluidModels();

    for (unsigned int modelIdx = 0; modelIdx < numFluidModels; modelIdx++)
    {
        std::vector<unsigned int> &cellSecondBorderLevels = m_cellSecondBorderLevels[modelIdx];
        std::fill_n(cellSecondBorderLevels.data(), cellSecondBorderLevels.size(), UINT32_MAX);

        std::vector<unsigned int>& particleSecondBorderLevels = m_particleSecondBorderLevels[modelIdx];
        std::fill_n(particleSecondBorderLevels.data(), particleSecondBorderLevels.size(), UINT32_MAX);
    }

    for (unsigned int modelIdx = 0; modelIdx < numFluidModels; modelIdx++)
    {
        #pragma omp parallel default(shared)
        {
            const std::vector<unsigned int>& cellRegionLevels = m_cellRegionLevels[modelIdx];
            const std::vector<unsigned int>& cellBorderLevels = m_cellBorderLevels[modelIdx];
            std::vector<unsigned int>& cellSecondBorderLevels = m_cellSecondBorderLevels[modelIdx];
            const std::vector<CellParticlePairIndices>& cellParticlePairIndices = m_cellPairIndices[modelIdx];
            std::vector<unsigned int>& particleSecondBorderLevels = m_particleSecondBorderLevels[modelIdx];
            const std::vector<CellParticlePair>& cellParticlePairs = m_cellParticlePairs[modelIdx];

            const int cellCount = (int)cellBorderLevels.size();
            const int zStep = m_resolution.x() * m_resolution.y();

            #pragma omp for schedule(static)
            for (int cellIdx = 0; cellIdx < cellCount; cellIdx++)
            {
                const CellParticlePairIndices& pairIndices = cellParticlePairIndices[cellIdx];

                const unsigned int cellLevel = cellRegionLevels[cellIdx];

                // Skip empty cells and level 0 cells
                if (pairIndices.pairIndexBegin == pairIndices.pairIndexEnd || cellLevel == 0)
                    continue;

                const unsigned int cellBorderLevel = cellBorderLevels[cellIdx];
                unsigned int& cellSecondBorderLevel = cellSecondBorderLevels[cellIdx];

                const int cellPosZ = cellIdx / zStep;
                const int cellPosY = (cellIdx - cellPosZ * zStep) / m_resolution.x();
                const int cellPosX = cellIdx % m_resolution.x();

                const int offsetStartX = std::max(-1, -cellPosX);
                const int offsetStartY = std::max(-1, -cellPosY);
                const int offsetStartZ = std::max(-1, -cellPosZ);

                const int offsetEndX = std::min(2, m_resolution.x() - cellPosX);
                const int offsetEndY = std::min(2, m_resolution.y() - cellPosY);
                const int offsetEndZ = std::min(2, m_resolution.z() - cellPosZ);

                // Check all neighbors in a 3x3x3 grid around the cell
                for (int offsetZ = offsetStartZ; offsetZ < offsetEndZ; offsetZ++)
                {
                    for (int offsetY = offsetStartY; offsetY < offsetEndY; offsetY++)
                    {
                        for (int offsetX = offsetStartX; offsetX < offsetEndX; offsetX++)
                        {
                            const unsigned int neighborIdx = cellIdx + offsetX + offsetY * m_resolution.x() + offsetZ * zStep;
                            const CellParticlePairIndices& neighborCellPairIndices = cellParticlePairIndices[neighborIdx];
                            const unsigned int neighborBorderLevel = cellBorderLevels[neighborIdx];

                            if (neighborBorderLevel < cellLevel)
                            {
                                cellSecondBorderLevel = std::min(cellSecondBorderLevel, neighborBorderLevel);
                            }
                        }
                    }
                }

                if (cellSecondBorderLevel != UINT32_MAX)
                {
                    unsigned int pairIdx = pairIndices.pairIndexBegin;
                    const unsigned int pairIdxEnd = pairIndices.pairIndexEnd;

                    while (pairIdx < pairIdxEnd)
                    {
                        const unsigned int particleIdx = cellParticlePairs[pairIdx++].particleIndex;
                        particleSecondBorderLevels[particleIdx] = cellSecondBorderLevel;
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
        std::vector<unsigned int>& secondBorderParticleIndices = m_borderParticleIndices[modelIdx];
        borderParticleIndices.clear();
        secondBorderParticleIndices.clear();

        const std::vector<unsigned int>& borderLevels = m_cellBorderLevels[modelIdx];
        const std::vector<unsigned int>& secondBorderLevels = m_cellSecondBorderLevels[modelIdx];
        const std::vector<CellParticlePairIndices>& cellParticlePairIndices = m_cellPairIndices[modelIdx];
        const std::vector<CellParticlePair>& cellParticlePairs = m_cellParticlePairs[modelIdx];

        const unsigned int cellCount = borderLevels.size();

        for (unsigned int level = 0; level < REGION_LEVELS_COUNT - 1; level++)
        {
            const unsigned int indicesInitialSize = borderParticleIndices.size();

            // Write border particle indices
            for (unsigned int cellIdx = 0; cellIdx < cellCount; cellIdx++)
            {
                if (borderLevels[cellIdx] == level)
                {
                    const CellParticlePairIndices& pairIndices = cellParticlePairIndices[cellIdx];

                    unsigned int pairIdx = pairIndices.pairIndexBegin;
                    const unsigned int pairIdxEnd = pairIndices.pairIndexEnd;

                    while (pairIdx < pairIdxEnd)
                    {
                        const unsigned int particleIdx = cellParticlePairs[pairIdx++].particleIndex;
                        borderParticleIndices.push_back(particleIdx);
                    }
                }
            }

            m_levelBorderParticleCounts[level][modelIdx] = borderParticleIndices.size() - indicesInitialSize;

            // Write border neighbor particle indices
            for (unsigned int cellIdx = 0; cellIdx < cellCount; cellIdx++)
            {
                if (secondBorderLevels[cellIdx] == level)
                {
                    const CellParticlePairIndices& pairIndices = cellParticlePairIndices[cellIdx];

                    unsigned int pairIdx = pairIndices.pairIndexBegin;
                    const unsigned int pairIdxEnd = pairIndices.pairIndexEnd;

                    while (pairIdx < pairIdxEnd)
                    {
                        const unsigned int particleIdx = cellParticlePairs[pairIdx++].particleIndex;
                        secondBorderParticleIndices.push_back(particleIdx);
                    }
                }
            }

            m_nSecondBorderUnionParticles[level][modelIdx] = secondBorderParticleIndices.size();
        }
    }
}
