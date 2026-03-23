#include <cmath>

#include "Scene.h"

#include "Numeric/PCGSolver.h"

namespace asa
{
    namespace
    {
        template<typename T> T interpolate(const Array2<T>& array, const Vector2 fractionalIndex)
        {
            auto i0 = std::floor(fractionalIndex.x);
            auto j0 = std::floor(fractionalIndex.y);
            auto i1 = i0 + 1;
            auto j1 = j0 + 1;

            const auto tx = fractionalIndex.x - i0;
            const auto ty = fractionalIndex.y - j0;

            i0 = clamp(i0, 0, array.getSize().x - 1);
            i1 = clamp(i1, 0, array.getSize().x - 1);
            j0 = clamp(j0, 0, array.getSize().y - 1);
            j1 = clamp(j1, 0, array.getSize().y - 1);

            return bilerp<T>(
                array.getValue(i0, j0),
                array.getValue(i1, j0),
                array.getValue(i0, j1),
                array.getValue(i1, j1),
                tx,
                ty
            );
        }

        SparseMatrix<float> A;
        void buildMatrixA(const uint sizeX, const uint sizeY, const float cellSizeX2, const float cellSizeY2)
        {
            if (A.n == sizeX * sizeY) return;

            A.clear();
            A.resize(sizeX * sizeY);

            const auto horizontalCoef = 1/cellSizeX2;
            const auto verticalCoef = 1/cellSizeY2;
            for (int i = 0; i < sizeX; ++i)
            {
                for (int j = 0; j < sizeY; ++j)
                {
                    const auto idx = i + j * sizeX;

                    A.set_element(idx, idx, 0.0f);

                    if (i + 1 < sizeX)
                    {
                        A.set_element(idx, idx + 1, - horizontalCoef);
                        A.add_to_element(idx, idx, horizontalCoef);
                    }

                    if (i - 1 >= 0)
                    {
                        A.set_element(idx, idx - 1, - horizontalCoef);
                        A.add_to_element(idx, idx, horizontalCoef);
                    }

                    if (j + 1 < sizeY)
                    {
                        A.set_element(idx, idx + sizeX, - verticalCoef);
                        A.add_to_element(idx, idx, verticalCoef);
                    }

                    if (j - 1 >= 0)
                    {
                        A.set_element(idx, idx - sizeX, - verticalCoef);
                        A.add_to_element(idx, idx, verticalCoef);
                    }
                }
            }

            A.set_element(0, 0, 1.0f);
            if (1 < sizeX * sizeY) A.set_element(0, 1, 0.0f);
            if (sizeX < sizeX * sizeY) A.set_element(0, sizeX, 0.0f);
        }
    }

    void Fluid2::fluidEmission()
    {
        if (Scene::testcase >= Scene::SMOKE)
        {
            for (int i = 0; i < grid.getSize().x; ++i)
            {
                for (int j = 0; j < grid.getSize().y; ++j)
                {
                    auto pos = grid.getCellPos(Index2(i, j));
                    
                    if (pos.x >= -0.1f && pos.x <= 0.1f && pos.y >= -1.9f && pos.y <= -1.75f)
                        inkRGB.setValue(i, j, Vector3(1.0f, 1.0f, 0.0f));
                    else if (pos.x >= -0.2f && pos.x < -0.1f && pos.y >= -1.9f && pos.y <= -1.75f)
                        inkRGB.setValue(i, j, Vector3(1.0f, 0.0f, 1.0f));
                    else if (pos.x > 0.1f && pos.x <= 0.2f && pos.y >= -1.9f && pos.y <= -1.75f)
                        inkRGB.setValue(i, j, Vector3(0.0f, 1.0f, 1.0f));
                }
            }
        }
    }

    void Fluid2::fluidAdvection(const float dt)
    {

        {
            // Ink advection HERE
            const auto oldInkRGB = inkRGB;

            for (int i = 0; i < grid.getSize().x; ++i)
            {
                for (int j = 0; j < grid.getSize().y; ++j)
                {
                    auto pos = grid.getCellPos(Index2(i, j));
                    auto vel = Vector2{
                        interpolate(velocityX, grid.getFaceIndexX(pos)),
                        interpolate(velocityY, grid.getFaceIndexY(pos))
                    };
                    auto prevPos = pos - dt * vel;
                    inkRGB.setValue(Index2(i, j), interpolate(oldInkRGB, grid.getCellIndex(prevPos)));
                }
            }
        }

        {
            // Velocity advection term HERE
            const auto oldU = velocityX;
            const auto oldV = velocityY;

            for (int i = 0; i < grid.getSizeFacesX().x; ++i)
            {
                for (int j = 0; j < grid.getSizeFacesX().y; ++j)
                {
                    auto pos = grid.getFacePosX(Index2(i, j));
                    auto vel = Vector2{
                        oldU.getValue(i, j),
                        interpolate(oldV, grid.getFaceIndexY(pos))
                    };
                    auto prevPos = pos - dt * vel;
                    velocityX.setValue(Index2(i, j), interpolate(oldU, grid.getFaceIndexX(prevPos)));
                }
            }

            for (int i = 0; i < grid.getSizeFacesY().x; ++i)
            {
                for (int j = 0; j < grid.getSizeFacesY().y; ++j)
                {
                    auto pos = grid.getFacePosY(Index2(i, j));
                    auto vel = Vector2{
                        interpolate(oldU, grid.getFaceIndexX(pos)),
                        oldV.getValue(i, j)
                    };
                    auto prevPos = pos - dt * vel;
                    velocityY.setValue(Index2(i, j), interpolate(oldV, grid.getFaceIndexY(prevPos)));
                }
            }
        }
    }

    void Fluid2::fluidVolumeForces(const float dt)
    {
        if (Scene::testcase >= Scene::SMOKE)
        {
            // External forces term HERE

            // Gravity
            for (int i = 0; i < grid.getSizeFacesY().x; ++i)
                for (int j = 0; j < grid.getSizeFacesY().y; ++j)
                    velocityY.setValue(Index2(i, j), velocityY.getValue(i, j) + dt * Scene::kGravity);

            // Emitters
            for (int i = 0; i < grid.getSizeFacesY().x; ++i)
            {
                for (int j = 0; j < grid.getSizeFacesY().y; ++j)
                {
                    auto pos = grid.getFacePosY(Index2(i, j));
                    if (pos.x >= -0.2f && pos.x <= 0.2f && pos.y >= -1.9f && pos.y <= -1.75f)
                    {
                        velocityY.setValue(Index2(i, j), 8.0f);
                    }
                }
            }
            
            for (int i = 0; i < grid.getSizeFacesX().x; ++i)
            {
                for (int j = 0; j < grid.getSizeFacesX().y; ++j)
                {
                    auto pos = grid.getFacePosX(Index2(i, j));
                    if (pos.x >= -0.2f && pos.x <= 0.2f && pos.y >= -1.9f && pos.y <= -1.75f)
                    {
                        velocityX.setValue(Index2(i, j), 0.0f);
                    }
                }
            }
        }
    }

    void Fluid2::fluidViscosity(const float dt)
    {
        if (Scene::testcase >= Scene::SMOKE)
        {
            // Viscosity term HERE
            const auto oldU = velocityX;
            const auto oldV = velocityY;
            const auto [cellSizeX, cellSizeY] = grid.getDx();
            const auto cellSizeX2 = cellSizeX * cellSizeX;
            const auto cellSizeY2 = cellSizeY * cellSizeY;
            const auto invCellSizeX2 = 1.0f / cellSizeX2;
            const auto invCellSizeY2 = 1.0f / cellSizeY2;

            for (int i = 1; i < grid.getSizeFacesX().x - 1; ++i)
            {
                for (int j = 1; j < grid.getSizeFacesX().y - 1; ++j)
                {
                    const auto u = Index2(i, j);
                    const auto u_i_j = oldU.getValue(u);
                    const auto u_ip1_j = oldU.getValue(Index2(i + 1, j));
                    const auto u_im1_j = oldU.getValue(Index2(i - 1, j));
                    const auto u_i_jp1 = oldU.getValue(Index2(i, j + 1));
                    const auto u_i_jm1 = oldU.getValue(Index2(i, j - 1));

                    auto laplacianHorizontal = u_ip1_j - 2 * u_i_j + u_im1_j;
                    laplacianHorizontal *= invCellSizeX2;
                    auto laplacianVertical = u_i_jp1 - 2 * u_i_j + u_i_jm1;
                    laplacianVertical *= invCellSizeY2;
                    const auto laplacian = laplacianHorizontal + laplacianVertical;

                    velocityX.setValue(u,  u_i_j + dt * Scene::kViscosity * laplacian);
                }
            }

            for (int i = 1; i < grid.getSizeFacesY().x - 1; ++i)
            {
                for (int j = 1; j < grid.getSizeFacesY().y - 1; ++j)
                {
                    const auto v = Index2(i, j);
                    const auto v_i_j = oldV.getValue(v);
                    const auto v_ip1_j = oldV.getValue(Index2(i + 1, j));
                    const auto v_im1_j = oldV.getValue(Index2(i - 1, j));
                    const auto v_i_jp1 = oldV.getValue(Index2(i, j + 1));
                    const auto v_i_jm1 = oldV.getValue(Index2(i, j - 1));

                    auto laplacianHorizontal = v_ip1_j - 2 * v_i_j + v_im1_j;
                    laplacianHorizontal /= cellSizeX2;
                    auto laplacianVertical = v_i_jp1 - 2 * v_i_j + v_i_jm1;
                    laplacianVertical /= cellSizeY2;
                    const auto laplacian = laplacianHorizontal + laplacianVertical;

                    velocityY.setValue(v,  v_i_j + dt * Scene::kViscosity * laplacian);
                }
            }
        }
    }

    void Fluid2::fluidPressureProjection(const float dt)
    {
        if (Scene::testcase >= Scene::SMOKE)
        {
            // Incompressibility / Pressure term HERE
            const auto [sizeX, sizeY] = grid.getSize();
            const auto [cellSizeX, cellSizeY] = grid.getDx();
            const auto cellSizeX2 = cellSizeX * cellSizeX;
            const auto cellSizeY2 = cellSizeY * cellSizeY;

            // Boundary conditions
            for (int j = 0; j < sizeY; ++j)
            {
                velocityX.setValue(Index2(0, j), 0.0f);
                velocityX.setValue(Index2(sizeX, j), 0.0f);
            }

            for (int i = 0; i < sizeX; ++i)
            {
                velocityY.setValue(Index2(i, 0), 0.0f);
                velocityY.setValue(Index2(i, sizeY), 0.0f);
            }

            // Build vector b
            static std::vector<float> b;
            b.assign(sizeX * sizeY, 0.0f);

            for (int i = 0; i < sizeX; ++i)
            {
                for (int j = 0; j < sizeY; ++j)
                {
                    const auto row = i + j * sizeX;

                    const auto left = velocityX.getValue(i, j);
                    const auto right = velocityX.getValue(i + 1, j);
                    const auto bottom = velocityY.getValue(i, j);
                    const auto top = velocityY.getValue(i, j + 1);

                    const auto discreteDivergence = (right - left)/cellSizeX + (top - bottom)/cellSizeY;

                    b[row] = - discreteDivergence * Scene::kDensity / dt;
                }
            }
            b[0] = 0.0f;

            // Build matrix A
            buildMatrixA(sizeX, sizeY, cellSizeX2, cellSizeY2);

            // Solve system
            static std::vector<float> x;
            x.assign(sizeX * sizeY, 0.0f);
            float residual;
            int its;

            auto solver = PCGSolver<float>();
            solver.solve(A, b, x, residual, its);

            // Copy into pressure matrix
            for (int i = 0; i < sizeX; ++i)
                for (int j = 0; j < sizeY; ++j)
                    pressure.setValue(i, j, x[i + j * sizeX]);

            // Apply pressure gradients
            for (int i = 1; i < grid.getSizeFacesX().x - 1; ++i)
            {
                for (int j = 0; j < grid.getSizeFacesX().y; ++j)
                {
                    const auto leftPressure = pressure.getValue(i - 1, j);
                    const auto rightPressure = pressure.getValue(i, j);

                    const auto pressureGradient = (rightPressure - leftPressure)/cellSizeX;
                    velocityX.setValue(Index2(i, j), velocityX.getValue(i, j) - pressureGradient *  dt / Scene::kDensity);
                }
            }

            for (int i = 0; i < grid.getSizeFacesY().x; ++i)
            {
                for (int j = 1; j < grid.getSizeFacesY().y - 1; ++j)
                {
                    const auto bottomPressure = pressure.getValue(i, j - 1);
                    const auto topPressure = pressure.getValue(i, j);

                    const auto pressureGradient = (topPressure - bottomPressure)/cellSizeY;
                    velocityY.setValue(Index2(i, j), velocityY.getValue(i, j) - pressureGradient *  dt / Scene::kDensity);
                }
            }


        }
    }
}
