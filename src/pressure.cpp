#include "pressure.h"
#include <Eigen/IterativeLinearSolvers>

void pressure(
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
    triplets.reserve(nx * ny * nz * 12 - 4 * (nx * ny + nx * nz + ny * nz));
    
    Eigen::VectorXd B;
    B.resize(6);
    B << -1.0 / grid.h, 1.0 / grid.h, -1.0 / grid.h, 1.0 / grid.h, -1.0 / grid.h, 1.0 / grid.h;
    
    Eigen::VectorXd q;
    q.resize(6);
    
    Eigen::VectorXd d;
    d.resize(nx * ny * nz);
    d.setZero();
    
    for (int i = 0; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            for (int k = 0; k < nz; k ++) {
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
                d(i + j * nx + k * nx * ny) = rho/dt * B.dot(q);
                
                //partial derivative
                //x derivative Pijk - P(i-1)jk
                if ((i - 1) >= 0){
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), ((i - 1) + j * nx + k * nx * ny), (-1.0 / grid.h) * (-1.0 / grid.h)));
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + j * nx + k * nx * ny), (1.0 / grid.h) * (-1.0 / grid.h)));
                }
                
                //x derivative P(i+1)jk - Pijk
                if ((i + 1) < nx) {
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + j * nx + k * nx * ny), (-1.0 / grid.h) * (1.0 / grid.h)));
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), ((i + 1) + j * nx + k * nx * ny), (1.0 / grid.h) * (1.0 / grid.h)));
                }
                
                //y derivative Pijk - Pi(j-1)k
                if ((j - 1) >= 0){
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + (j - 1) * nx + k * nx * ny), (-1.0 / grid.h) * (-1.0 / grid.h)));
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + j * nx + k * nx * ny), (1.0 / grid.h) * (-1.0 / grid.h)));
                }
                
                //y derivative Pi(j+1)k - Pijk
                if ((j + 1) < ny) {
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + j * nx + k * nx * ny), (-1.0 / grid.h) * (1.0 / grid.h)));
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + (j + 1) * nx + k * nx * ny), (1.0 / grid.h) * (1.0 / grid.h)));
                }
                
                //z derivative Pijk - Pij(k-1)
                if ((k - 1) >= 0){
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + j * nx + (k - 1) * nx * ny), (-1.0 / grid.h) * (-1.0 / grid.h)));
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + j * nx + k * nx * ny), (1.0 / grid.h) * (-1.0 / grid.h)));
                }
                
                //z derivative Pij(k+1) - Pijk
                if ((k + 1) < nz){
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + j * nx + k * nx * ny), (-1.0 / grid.h) * (1.0 / grid.h)));
                    triplets.push_back(Eigen::Triplet<double>((i + j * nx + k * nx * ny), (i + j * nx + (k + 1) * nx * ny), (1.0 / grid.h) * (1.0 / grid.h)));
                }
            }
        }
    }
    //construct the global D matrix
    Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>(nx * ny * nz, nx * ny * nz);
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    grid.P.resize(nx * ny * nz);
    grid.P.setZero();
    ConjugateGradient<Eigen::SparseMatrix<double>, Lower|Upper> cg;
    cg.compute(A);
    grid.P = cg.solve(d);
    
}
