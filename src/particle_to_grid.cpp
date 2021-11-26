#include "particle_to_grid.h"
#include "p_to_grid_interpolate.h"
#include <iostream>
void particle_to_grid(
    grid_data & grid,
    const Eigen::MatrixXd & Position,
    const Eigen::MatrixXd & Velocity,
    const Eigen::RowVector3d & corner,
    const int nx,
    const int ny,
    const int nz,
    const double h)
{
    //initialze constant;
    grid.nx = nx;
    grid.ny = ny;
    grid.nz = nz;
    grid.h = h;
    grid.corner = corner;
    Eigen::VectorXd Grid;
    //use trilinear interpolation to get X, Y, Z
    Eigen::RowVector3d temp_corner;
    //0 X
    temp_corner = corner;
    temp_corner(1) += 0.5 * h;
    temp_corner(2) += 0.5 * h;
    
    Grid.resize((nx + 1) * ny * nz);
    Grid.setZero();
    p_to_grid_interpolate(nx + 1, ny, nz, h, temp_corner, Position, Velocity.col(0), 0, Grid);
    grid.X = Grid;
    //1 Y
    temp_corner = corner;
    temp_corner(0) += 0.5 * h;
    temp_corner(2) += 0.5 * h;
    
    Grid.resize(nx * (ny + 1) * nz);
    Grid.setZero();
    p_to_grid_interpolate(nx, ny + 1, nz, h, temp_corner, Position, Velocity.col(1), 1, Grid);
    grid.Y = Grid;
    
    //2 Z
    temp_corner = corner;
    temp_corner(0) += 0.5 * h;
    temp_corner(1) += 0.5 * h;
    
    Grid.resize(nx * ny * (nz + 1));
    Grid.setZero();
    p_to_grid_interpolate(nx, ny, nz + 1, h, temp_corner, Position, Velocity.col(2), 2, Grid);
    grid.Z = Grid;
    
    grid.P.resize(nx * ny * nz);
    grid.P.setZero();
}
