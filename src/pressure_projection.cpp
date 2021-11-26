#include "pressure_projection.h"
#include <Eigen/IterativeLinearSolvers>

void pressure_projection(
    grid_data & grid,
    const double rho,
    const double dt)
{
    
    int nx, ny, nz;
    nx = grid.nx;
    ny = grid.ny;
    nz = grid.nz;
    
    //separate version
    // X
    for (int i = 1; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            for (int k = 0; k < nz; k ++) {
                grid.X(i + j * (nx + 1) + k * (nx + 1) * ny) = (-1.0 * dt / rho) * (grid.P(i + j * nx + k * nx * ny) - grid.P((i - 1) + j * nx + k * nx * ny)) / grid.h;
            }
        }
    }
    
    // Y
    for (int i = 0; i < nx; i ++) {
        for (int j = 1; j < ny; j ++) {
            for (int k = 0; k < nz; k ++) {
                grid.Y(i + j * nx + k * nx * (ny + 1)) = (-1.0 * dt / rho) * (grid.P(i + j * nx + k * nx * ny) - grid.P(i + (j - 1) * nx + k * nx * ny)) / grid.h;
            }
        }
    }
    
    // Z
    for (int i = 0; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            for (int k = 1; k < nz; k ++) {
                grid.Z(i + j * nx + k * nx * ny) = (-1.0 * dt / rho) * (grid.P(i + j * nx + k * nx * ny) - grid.P(i + j * nx + (k - 1) * nx * ny)) / grid.h;

            }
        }
    }
    
    //compact version
    //for (int i = 0; i < nx; i ++) {
    //    for (int j = 0; j < ny; j ++) {
    //        for (int k = 0; k < nz; k ++) {
    //            if (i > 0) {
    //                grid.X(i + j * (nx + 1) + k * (nx + 1) * ny) = (-1.0 * dt / rho) * (grid.P(i + j * nx + k * nx * ny) - grid.P((i - 1) + j * nx + k * nx * ny)) / grid.h;
    //            }
    //            if (j > 0) {
    //                grid.Y(i + j * nx + k * nx * (ny + 1)) = (-1.0 * dt / rho) * (grid.P(i + j * nx + k * nx * ny) - grid.P(i + (j - 1) * nx + k * nx * ny)) / grid.h;
    //            }
    //            if (k > 0) {
    //                grid.Z(i + j * nx + k * nx * ny) = (-1.0 * dt / rho) * (grid.P(i + j * nx + k * nx * ny) - grid.P(i + j * nx + (k - 1) * nx * ny)) / grid.h;
    //            }
    //        }
    //    }
    //}
}
