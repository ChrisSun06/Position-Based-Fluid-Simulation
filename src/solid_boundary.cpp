#include "solid_boundary.h"
void solid_boundary(
    grid_data & grid)
{
    //X i
    for (int j = 0; j < grid.ny; j ++) {
        for (int k = 0; k < grid.nz; k ++) {
            grid.X(      0 + j * (grid.nx + 1) + k * (grid.nx + 1) * grid.ny) = -1.0 * grid.X(      0 + j * (grid.nx + 1) + k * (grid.nx + 1) * grid.ny);
            grid.X(grid.nx + j * (grid.nx + 1) + k * (grid.nx + 1) * grid.ny) = -1.0 * grid.X(grid.nx + j * (grid.nx + 1) + k * (grid.nx + 1) * grid.ny);
        }
    }
    //Y j
    for (int i = 0; i < grid.nx; i ++) {
        for (int k = 0; k < grid.nz; k ++) {
            grid.Y(i +       0 * grid.nx + k * grid.nx * (grid.ny + 1)) = -1.0 * grid.Y(i +       0 * grid.nx + k * grid.nx * (grid.ny + 1));
            grid.Y(i + grid.ny * grid.nx + k * grid.nx * (grid.ny + 1)) = -1.0 * grid.Y(i + grid.ny * grid.nx + k * grid.nx * (grid.ny + 1));
        }
    }
    //Z k
    for (int i = 0; i < grid.nx; i ++) {
        for (int j = 0; j < grid.ny; j ++) {
            grid.Z(i + j * grid.nx +       0 * grid.nx * grid.ny) = -1.0 * grid.Z(i + j * grid.nx +       0 * grid.nx * grid.ny);
            grid.Z(i + j * grid.nx + grid.nz * grid.nx * grid.ny) = -1.0 * grid.Z(i + j * grid.nx + grid.nz * grid.nx * grid.ny);
        }
    }
    
}
