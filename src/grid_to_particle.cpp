#include "grid_to_particle.h"
#include "grid_to_p_interpolate.h"
void grid_to_particle(
    const grid_data & grid,
    const Eigen::MatrixXd & Position,
    Eigen::MatrixXd & Velocity)
{
    //initialze constant;
    
    int nx, ny, nz;
    double h;
    Eigen::RowVector3d corner;
    nx = grid.nx;
    ny = grid.ny;
    nz = grid.nz;
    h = grid.h;
    corner = grid.corner;
    
    //use trilinear interpolation to get X, Y, Z
    Eigen::RowVector3d temp_corner;
    Eigen::VectorXd V;
    
    //0 X
    temp_corner = corner;
    temp_corner(1) += 0.5 * h;
    temp_corner(2) += 0.5 * h;
    V = Velocity.col(0);
    grid_to_p_interpolate(nx + 1, ny, nz, h, temp_corner, Position, V, 0, grid.X);
    Velocity.col(0) = V;
    
    //1 Y
    temp_corner = corner;
    temp_corner(0) += 0.5 * h;
    temp_corner(2) += 0.5 * h;
    V = Velocity.col(1);
    grid_to_p_interpolate(nx, ny + 1, nz, h, temp_corner, Position, V, 1, grid.Y);
    Velocity.col(1) = V;
    
    //2 Z
    temp_corner = corner;
    temp_corner(0) += 0.5 * h;
    temp_corner(1) += 0.5 * h;
    V = Velocity.col(2);
    grid_to_p_interpolate(nx, ny, nz + 1, h, temp_corner, Position, V, 2, grid.Z);
    Velocity.col(2) = V;
    
}
