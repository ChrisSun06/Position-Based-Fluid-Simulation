#include "smooth_air_fluid.h"
#include <cmath>

double Gaussian_3d(
  const double position,
  const double variance)
{
    return exp(-1.0 * position / (2.0 * variance));
}

void smooth_air_fluid(
  grid_data & grid)
{
    int nx, ny, nz;
    nx = grid.nx;
    ny = grid.ny;
    nz = grid.nz;
    
    //initialize all cells to air cells
    grid.air_fluid.resize(nx * ny * nz);
    //grid.air_fluid.setOnes();
    //grid.air_fluid = (grid.Raw_air_fluid * 1.0).eval();
    
    grid.fluid_index.resize(nx * ny * nz);
    grid.fluid_index.setOnes();
    grid.fluid_index = (grid.fluid_index * -1).eval();
    
    double weight;
    double value;
    
    double variance = 0.5;
    double d0 = Gaussian_3d(0.0, variance);
    double d1 = Gaussian_3d(1.0, variance);
    double d2 = Gaussian_3d(2.0, variance);
    double d3 = Gaussian_3d(3.0, variance);
    
    int counter = 0;
    grid.num_fluid_cells = 0;
    for (int i = 0; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            for (int k = 0; k < nz; k ++) {
                weight = 0.0;
                value = 0.0;
                
                //the centering cell;
                weight += d0;
                value += d0 * grid.Raw_air_fluid(i + j * nx + k * nx * ny);
                
                //check cells with distance 1, there are 6 of them
                if ((i - 1) >= 0){
                    weight += d1;
                    value += d1 * grid.Raw_air_fluid((i - 1) + j * nx + k * nx * ny);
                }
                if ((i + 1) < nx) {
                    weight += d1;
                    value += d1 * grid.Raw_air_fluid((i + 1) + j * nx + k * nx * ny);
                }
                if ((j - 1) >= 0){
                    weight += d1;
                    value += d1 * grid.Raw_air_fluid(i + (j - 1) * nx + k * nx * ny);
                }
                if ((j + 1) < ny) {
                    weight += d1;
                    value += d1 * grid.Raw_air_fluid(i + (j + 1) * nx + k * nx * ny);
                }
                if ((k - 1) >= 0){
                    weight += d1;
                    value += d1 * grid.Raw_air_fluid(i + j * nx + (k - 1) * nx * ny);
                }
                if ((k + 1) < nz) {
                    weight += d1;
                    value += d1 * grid.Raw_air_fluid(i + j * nx + (k + 1) * nx * ny);
                }
                
                //check cells with distance 2, there are 12 of them
                if ((i - 1) >= 0 && (j - 1) >= 0){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid((i - 1) + (j - 1) * nx + k * nx * ny);
                }
                if ((i - 1) >= 0 && (j + 1) < ny){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid((i - 1) + (j + 1) * nx + k * nx * ny);
                }
                if ((i + 1) < nx && (j - 1) >= 0){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid((i + 1) + (j - 1) * nx + k * nx * ny);
                }
                if ((i + 1) < nx && (j + 1) < ny){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid((i + 1) + (j + 1) * nx + k * nx * ny);
                }
                
                if ((i - 1) >= 0 && (k - 1) >= 0){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid((i - 1) + j * nx + (k - 1) * nx * ny);
                }
                if ((i - 1) >= 0 && (k + 1) < nz){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid((i - 1) + j * nx + (k + 1) * nx * ny);
                }
                if ((i + 1) < nx && (k - 1) >= 0){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid((i + 1) + j * nx + (k - 1) * nx * ny);
                }
                if ((i + 1) < nx && (k + 1) < nz){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid((i + 1) + j * nx + (k + 1) * nx * ny);
                }
                
                if ((j - 1) >= 0 && (k - 1) >= 0){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid(i + (j - 1) * nx + (k - 1) * nx * ny);
                }
                if ((j - 1) >= 0 && (k + 1) < nz){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid(i + (j - 1) * nx + (k + 1) * nx * ny);
                }
                if ((j + 1) < ny && (k - 1) >= 0){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid(i + (j + 1) * nx + (k - 1) * nx * ny);
                }
                if ((j + 1) < ny && (k + 1) < nz){
                    weight += d2;
                    value += d2 * grid.Raw_air_fluid(i + (j + 1) * nx + (k + 1) * nx * ny);
                }
                
                
                //check cells with distance 3, there are 8 of them
                if ((i - 1) >= 0 && (j - 1) >= 0 && (k - 1) >= 0){
                    weight += d3;
                    value += d3 * grid.Raw_air_fluid((i - 1) + (j - 1) * nx + (k - 1) * nx * ny);
                }
                if ((i - 1) >= 0 && (j - 1) >= 0 && (k + 1) < nz){
                    weight += d3;
                    value += d3 * grid.Raw_air_fluid((i - 1) + (j - 1) * nx + (k + 1) * nx * ny);
                }
                if ((i - 1) >= 0 && (j + 1) < ny && (k - 1) >= 0){
                    weight += d3;
                    value += d3 * grid.Raw_air_fluid((i - 1) + (j + 1) * nx + (k - 1) * nx * ny);
                }
                if ((i - 1) >= 0 && (j + 1) < ny && (k + 1) < nz){
                    weight += d3;
                    value += d3 * grid.Raw_air_fluid((i - 1) + (j + 1) * nx + (k + 1) * nx * ny);
                }
                if ((i + 1) < nx && (j - 1) >= 0 && (k - 1) >= 0){
                    weight += d3;
                    value += d3 * grid.Raw_air_fluid((i + 1) + (j - 1) * nx + (k - 1) * nx * ny);
                }
                if ((i + 1) < nx && (j - 1) >= 0 && (k + 1) < nz){
                    weight += d3;
                    value += d3 * grid.Raw_air_fluid((i + 1) + (j - 1) * nx + (k + 1) * nx * ny);
                }
                if ((i + 1) < nx && (j + 1) < ny && (k - 1) >= 0){
                    weight += d3;
                    value += d3 * grid.Raw_air_fluid((i + 1) + (j + 1) * nx + (k - 1) * nx * ny);
                }
                if ((i + 1) < nx && (j + 1) < ny && (k + 1) < nz){
                    weight += d3;
                    value += d3 * grid.Raw_air_fluid((i + 1) + (j + 1) * nx + (k + 1) * nx * ny);
                }
                
                grid.air_fluid(i + j * nx + k * nx * ny) = value / weight;
                //if value is less than 0, it is a fluid cell;
                if (value < 0.0) {
                    grid.num_fluid_cells = grid.num_fluid_cells + 1;
                    grid.fluid_index(i + j * nx + k * nx * ny) = counter;
                    counter += 1;
                }
            }
        }
    }
}

