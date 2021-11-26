#include "pressure_with_free_surface.h"
#include "weight_calculator.h"
#include <Eigen/IterativeLinearSolvers>

void pressure_with_free_surface(
    grid_data & grid,
    const double rho,
    const double dt)
{
    using namespace Eigen;
    
    int nx, ny, nz;
    nx = grid.nx;
    ny = grid.ny;
    nz = grid.nz;
    
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(grid.num_fluid_cells * 12);
    
    Eigen::VectorXd B;
    B.resize(6);
    B << -1.0 / grid.h, 1.0 / grid.h, -1.0 / grid.h, 1.0 / grid.h, -1.0 / grid.h, 1.0 / grid.h;
    
    Eigen::VectorXd q;
    q.resize(6);
    
    Eigen::VectorXd d;
    d.resize(grid.num_fluid_cells);
    d.setZero();
    
    int current_index, neighbouring_index;
    double low;
    double high;
    double alpha;
    int counter = 0;
    for (int i = 0; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            for (int k = 0; k < nz; k ++) {
                //only consider fluid cell;
                current_index = grid.fluid_index(i + j * nx + k * nx * ny);
                if (current_index > -1){
                    //construct q
                    q << grid.X(i + j * (nx + 1) + k * (nx + 1) * ny),
                         grid.X((i + 1) + j * (nx + 1) + k * (nx + 1) * ny),
                         grid.Y(i + j * nx + k * nx * (ny + 1)),
                         grid.Y(i + (j + 1) * nx + k * nx * (ny + 1)),
                         grid.Z(i + j * nx + k * nx * ny),
                         grid.Z(i + j * nx + (k + 1) * nx * ny);
                    
                    //solid boundary
                    if (i == 0) {
                        q(0) = 0.0;
                    }
                    if ((i + 1) == nx) {
                        q(1) = 0.0;
                    }
                
                    if (j == 0) {
                        q(2) = 0.0;
                    }
                    if ((j + 1) == ny) {
                        q(3) = 0.0;
                    }
                
                    if (k == 0) {
                        q(4) = 0.0;
                    }
                    if ((k + 1) == nz) {
                        q(5) = 0.0;
                    }
                
                    //construct d for global system;
                    d(counter) = rho/dt * B.dot(q);
                
                    //partial derivative
                    //x derivative Pijk - P(i-1)jk
                    //checking for solid boundary
                    if ((i - 1) >= 0){
                        neighbouring_index = grid.fluid_index((i - 1) + j * nx + k * nx * ny);
                        //check for free surface
                        if (neighbouring_index > -1) {
                            triplets.push_back(Eigen::Triplet<double>(counter, neighbouring_index, (-1.0 / grid.h) * (-1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (1.0 / grid.h) * (-1.0 / grid.h)));
                        } else {
                            //there is free surface
                            low = grid.air_fluid(i + j * nx + k * nx * ny);
                            high = grid.air_fluid((i - 1) + j * nx + k * nx * ny);
                            alpha = 1.0 - weight_calculator(low, high, 0.0);
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (1.0 / grid.h) * (-1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (alpha / (1.0 - alpha)) * (1.0 / grid.h) * (-1.0 / grid.h)));
                        }
                    } else {
                        //neighbouring cell is wall
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                    }
                
                    //x derivative P(i+1)jk - Pijk
                    //checking for solid boundary
                    if ((i + 1) < nx) {
                        neighbouring_index = grid.fluid_index((i + 1) + j * nx + k * nx * ny);
                        //check for free surface
                        if (neighbouring_index > -1) {
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (-1.0 / grid.h) * (1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, neighbouring_index, (1.0 / grid.h) * (1.0 / grid.h)));
                        } else {
                            //there is free surface
                            low = grid.air_fluid(i + j * nx + k * nx * ny);
                            high = grid.air_fluid((i + 1) + j * nx + k * nx * ny);
                            alpha = 1.0 - weight_calculator(low, high, 0.0);
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (-1.0 / grid.h) * (1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (alpha / (1.0 - alpha)) * (-1.0 / grid.h) * (1.0 / grid.h)));
                        }
                    } else {
                        //neighbouring cell is wall
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                    }
                
                    //y derivative Pijk - Pi(j-1)k
                    //checking for solid boundary
                    if ((j - 1) >= 0){
                        neighbouring_index = grid.fluid_index(i + (j - 1) * nx + k * nx * ny);
                        //check for free surface
                        if (neighbouring_index > -1) {
                            triplets.push_back(Eigen::Triplet<double>(counter, neighbouring_index, (-1.0 / grid.h) * (-1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (1.0 / grid.h) * (-1.0 / grid.h)));
                        } else {
                            //there is free surface
                            low = grid.air_fluid(i + j * nx + k * nx * ny);
                            high = grid.air_fluid(i + (j - 1) * nx + k * nx * ny);
                            alpha = 1.0 - weight_calculator(low, high, 0.0);
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (1.0 / grid.h) * (-1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (alpha / (1.0 - alpha)) * (1.0 / grid.h) * (-1.0 / grid.h)));
                        }
                    } else {
                        //neighbouring cell is wall
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                    }
                
                    //y derivative Pi(j+1)k - Pijk
                    //checking for solid boundary
                    if ((j + 1) < ny) {
                        neighbouring_index = grid.fluid_index(i + (j + 1) * nx + k * nx * ny);
                        //check for free surface
                        if (neighbouring_index > -1) {
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (-1.0 / grid.h) * (1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, neighbouring_index, (1.0 / grid.h) * (1.0 / grid.h)));
                        } else {
                            //there is free surface
                            low = grid.air_fluid(i + j * nx + k * nx * ny);
                            high = grid.air_fluid(i + (j + 1) * nx + k * nx * ny);
                            alpha = 1.0 - weight_calculator(low, high, 0.0);
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (-1.0 / grid.h) * (1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (alpha / (1.0 - alpha)) * (-1.0 / grid.h) * (1.0 / grid.h)));
                        }
                    } else {
                        //neighbouring cell is wall
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                    }
                
                    //z derivative Pijk - Pij(k-1)
                    //checking for solid boundary
                    if ((k - 1) >= 0){
                        neighbouring_index = grid.fluid_index(i + j * nx + (k - 1) * nx * ny);
                        //check for free surface
                        if (neighbouring_index > -1) {
                            triplets.push_back(Eigen::Triplet<double>(counter, neighbouring_index, (-1.0 / grid.h) * (-1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (1.0 / grid.h) * (-1.0 / grid.h)));
                        } else {
                            //there is free surface
                            low = grid.air_fluid(i + j * nx + k * nx * ny);
                            high = grid.air_fluid(i + j * nx + (k - 1) * nx * ny);
                            alpha = 1.0 - weight_calculator(low, high, 0.0);
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (1.0 / grid.h) * (-1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (alpha / (1.0 - alpha)) * (1.0 / grid.h) * (-1.0 / grid.h)));
                        }
                    } else {
                        //neighbouring cell is wall
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                    }
                
                    //z derivative Pij(k+1) - Pijk
                    //checking for solid boundary
                    if ((k + 1) < nz){
                        neighbouring_index = grid.fluid_index(i + j * nx + (k + 1) * nx * ny);
                        //check for free surface
                        if (neighbouring_index > -1) {
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (-1.0 / grid.h) * (1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, neighbouring_index, (1.0 / grid.h) * (1.0 / grid.h)));
                        } else {
                            //there is free surface
                            low = grid.air_fluid(i + j * nx + k * nx * ny);
                            high = grid.air_fluid(i + j * nx + (k + 1) * nx * ny);
                            alpha = 1.0 - weight_calculator(low, high, 0.0);
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (-1.0 / grid.h) * (1.0 / grid.h)));
                            triplets.push_back(Eigen::Triplet<double>(counter, current_index, (alpha / (1.0 - alpha)) * (-1.0 / grid.h) * (1.0 / grid.h)));
                        }
                    } else {
                        //neighbouring cell is wall
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                        triplets.push_back(Eigen::Triplet<double>(counter, current_index, 0.0));
                    }
                    
                    counter += 1;
                }
            }
        }
    }
    //construct the global D matrix
    Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>(grid.num_fluid_cells, grid.num_fluid_cells);
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    Eigen::VectorXd fluid_P;
    fluid_P.resize(grid.num_fluid_cells);
    fluid_P.setZero();
    
    ConjugateGradient<Eigen::SparseMatrix<double>, Lower|Upper> cg;
    cg.compute(A);
    fluid_P = cg.solve(d);
    
    
    grid.P.resize(nx * ny * nz);
    grid.P.setZero();
    for (int i = 0; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            for (int k = 0; k < nz; k ++) {
                current_index = grid.fluid_index(i + j * nx + k * nx * ny);
                if (current_index > -1){
                    grid.P(i + j * nx + k * nx * ny) = fluid_P(current_index);
                }
            }
        }
    }
}
