#include "raw_air_fluid.h"
#include <cmath>

void raw_air_fluid(
  const Eigen::MatrixXd & P,
  grid_data & grid)
{
    int i, j, k;
    
    int nx, ny, nz;
    nx = grid.nx;
    ny = grid.ny;
    nz = grid.nz;
    
    double h = grid.h;
    Eigen::RowVector3d corner = grid.corner;
    
    //initialize all cells to air cells
    grid.Raw_air_fluid.resize(nx * ny * nz);
    grid.Raw_air_fluid.setOnes();
    grid.Raw_air_fluid = (grid.Raw_air_fluid * 1.0).eval();
    
    
    for(int row = 0; row < P.rows(); row ++) {
        // find the cell in which the current particle stays.
        i = std::floor((P(row, 0) - corner(0)) / h);
        
        j = std::floor((P(row, 1) - corner(1)) / h);
        
        k = std::floor((P(row, 2) - corner(2)) / h);
        
        grid.Raw_air_fluid(i + j * nx + k * nx * ny) = -1.0;
        
    }

}

