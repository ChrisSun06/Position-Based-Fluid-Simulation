#ifndef GRID_DATA
#define GRID_DATA

#include <Eigen/Core>
#include <Eigen/Sparse>
struct grid_data
{
    int nx;
    int ny;
    int nz;
    
    //the size of one cell of the grid;
    double h;
    
    //starting corner of the regular grid, lower, left, back
    Eigen::RowVector3d corner;
    //stagger grid for Velocityx, size = (nx + 1) * ny * nz;
    //To access the element at i, j, k, we access the element with index i + j * (nx + 1) + k * (nx + 1) * ny
    Eigen::VectorXd X;
    //stagger grid for Velocityy, size = nx * (ny + 1) * nz;
    //To access the element at i, j, k, we access the element with index i + j * nx + k * nx * (ny + 1)
    Eigen::VectorXd Y;
    //stagger grid for Velocityz, size = nx * ny * (nz + 1);
    //To access the element at i, j, k, we access the element with index i + j * nx + k * nx * ny
    Eigen::VectorXd Z;
    //grid for pressure
    //To access the element at i, j, k, we access the element with index i + j * nx + k * nx * ny
    Eigen::VectorXd P;
    
    //Raw air_fluid_cell if a cell contains fluid particles, it is fluid, use -1 to represent. Using +1 to represent air cell, size nx * ny * nz;
    //To access the element at i, j, k, we access the element with index i + j * nx + k * nx * ny
    Eigen::VectorXd Raw_air_fluid;
    
    //smoothed version Raw_air_fluid, obtained by applying 3d convolution with gaussian kernel to Raw_air_fluid, size nx * ny * nz;
    //To access the element at i, j, k, we access the element with index i + j * nx + k * nx * ny
    Eigen::VectorXd air_fluid;
    
    int num_fluid_cells;
    //each fluid cell would have a unique index. air cells' index would be -1
    //To access the element at i, j, k, we access the element with index i + j * nx + k * nx * ny
    Eigen::VectorXi fluid_index;
};

#endif
