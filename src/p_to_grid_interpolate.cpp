#include "p_to_grid_interpolate.h"
#include "weight_calculator.h"
#include <cmath>

void p_to_grid_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  const Eigen::VectorXd & V,
  const int direction,
  Eigen::VectorXd & Grid)
{
    int x_low, x_high, y_low, y_high, z_low, z_high;
    double x_low_v, x_high_v, y_low_v, y_high_v, z_low_v, z_high_v;
    
    Grid.resize(nx * ny * nz);
    Grid.setZero();
    Grid = Eigen::VectorXd::Zero(nx * ny * nz);
    
    //Eigen::VectorXd Grid_weight;
    //Grid_weight.resize(nx * ny * nz);
    //Grid_weight.setZero();
    
    for(int row = 0; row < P.rows(); row ++) {
        // find the eight corners for interpolation.
        x_low = std::floor((P(row, 0) - corner(0)) / h);
        x_high = x_low + 1;
        x_low_v = corner(0) + h * x_low;
        x_high_v = x_low_v + h;
        
        y_low = std::floor((P(row, 1) - corner(1)) / h);
        y_high = y_low + 1;
        y_low_v = corner(1) + h * y_low;
        y_high_v = y_low_v + h;
        
        z_low = std::floor((P(row, 2) - corner(2)) / h);
        z_high = z_low + 1;
        z_low_v = corner(2) + h * z_low;
        z_high_v = z_low_v + h;
        
        x_high_v = weight_calculator(x_low_v, x_high_v, P(row, 0));
        x_low_v = 1 - x_high_v;
        
        y_high_v = weight_calculator(y_low_v, y_high_v, P(row, 1));
        y_low_v = 1 - y_high_v;
        
        z_high_v = weight_calculator(z_low_v, z_high_v, P(row, 2));
        z_low_v = 1 - z_high_v;

        //valid low index need to be greater than or equal to 0
        //valid high value should be lower than the number of node on corresponding axis
        //direction x, need to check direction y and z
        if (direction == 0) {
            if (z_low >= 0) {
                if (y_low >= 0){
                    Grid(x_low + y_low * nx + z_low * nx * ny) += x_low_v * y_low_v * z_low_v * V(row);
                    Grid(x_high + y_low * nx + z_low * nx * ny) += x_high_v * y_low_v * z_low_v * V(row);
                    
                    //Grid_weight(x_low + y_low * nx + z_low * nx * ny) += x_low_v * y_low_v * z_low_v;
                    //Grid_weight(x_high + y_low * nx + z_low * nx * ny) += x_high_v * y_low_v * z_low_v;
                }
                if (y_high < ny){
                    Grid(x_low + y_high * nx + z_low * nx * ny) += x_low_v * y_high_v * z_low_v * V(row);
                    Grid(x_high + y_high * nx + z_low * nx * ny) += x_high_v * y_high_v * z_low_v * V(row);
                    
                    //Grid_weight(x_low + y_high * nx + z_low * nx * ny) += x_low_v * y_high_v * z_low_v;
                    //Grid_weight(x_high + y_high * nx + z_low * nx * ny) += x_high_v * y_high_v * z_low_v;
                }
            }
            
            if (z_high < nz) {
                if (y_low >= 0){
                    Grid(x_low + y_low * nx + z_high * nx * ny) += x_low_v * y_low_v * z_high_v * V(row);
                    Grid(x_high + y_low * nx + z_high * nx * ny) += x_high_v * y_low_v * z_high_v * V(row);
                    
                    //Grid_weight(x_low + y_low * nx + z_high * nx * ny) += x_low_v * y_low_v * z_high_v;
                    //Grid_weight(x_high + y_low * nx + z_high * nx * ny) += x_high_v * y_low_v * z_high_v;
                }
                if (y_high < ny){
                    Grid(x_low + y_high * nx + z_high * nx * ny) += x_low_v * y_high_v * z_high_v * V(row);
                    Grid(x_high + y_high * nx + z_high * nx * ny) += x_high_v * y_high_v * z_high_v * V(row);
                    
                    //Grid_weight(x_low + y_high * nx + z_high * nx * ny) += x_low_v * y_high_v * z_high_v;
                    //Grid_weight(x_high + y_high * nx + z_high * nx * ny) += x_high_v * y_high_v * z_high_v;
                }
            }
        }
        
        //direction y, need to check direction x and z
        if (direction == 1) {
            if (z_low >= 0) {
                if (x_low >= 0){
                    Grid(x_low + y_low * nx + z_low * nx * ny) += x_low_v * y_low_v * z_low_v * V(row);
                    Grid(x_low + y_high * nx + z_low * nx * ny) += x_low_v * y_high_v * z_low_v * V(row);
                    
                    //Grid_weight(x_low + y_low * nx + z_low * nx * ny) += x_low_v * y_low_v * z_low_v;
                    //Grid_weight(x_low + y_high * nx + z_low * nx * ny) += x_low_v * y_high_v * z_low_v;
                }
                if (x_high < nx){
                    Grid(x_high + y_low * nx + z_low * nx * ny) += x_high_v * y_low_v * z_low_v * V(row);
                    Grid(x_high + y_high * nx + z_low * nx * ny) += x_high_v * y_high_v * z_low_v * V(row);
                    
                    //Grid_weight(x_high + y_low * nx + z_low * nx * ny) += x_high_v * y_low_v * z_low_v;
                    //Grid_weight(x_high + y_high * nx + z_low * nx * ny) += x_high_v * y_high_v * z_low_v;
                }
            }
            
            if (z_high < nz) {
                if (x_low >= 0){
                    Grid(x_low + y_low * nx + z_high * nx * ny) += x_low_v * y_low_v * z_high_v * V(row);
                    Grid(x_low + y_high * nx + z_high * nx * ny) += x_low_v * y_high_v * z_high_v * V(row);
                    
                    //Grid_weight(x_low + y_low * nx + z_high * nx * ny) += x_low_v * y_low_v * z_high_v;
                    //Grid_weight(x_low + y_high * nx + z_high * nx * ny) += x_low_v * y_high_v * z_high_v;
                }
                if (x_high < nx){
                    Grid(x_high + y_low * nx + z_high * nx * ny) += x_high_v * y_low_v * z_high_v * V(row);
                    Grid(x_high + y_high * nx + z_high * nx * ny) += x_high_v * y_high_v * z_high_v * V(row);
                    
                    //Grid_weight(x_high + y_low * nx + z_high * nx * ny) += x_high_v * y_low_v * z_high_v;
                    //Grid_weight(x_high + y_high * nx + z_high * nx * ny) += x_high_v * y_high_v * z_high_v;
                }
            }
        }

        //direction z, need to check direction x and y
        if (direction == 2) {
            if (x_low >= 0) {
                if (y_low >= 0){
                    Grid(x_low + y_low * nx + z_low * nx * ny) += x_low_v * y_low_v * z_low_v * V(row);
                    Grid(x_low + y_low * nx + z_high * nx * ny) += x_low_v * y_low_v * z_high_v * V(row);
                    
                    //Grid_weight(x_low + y_low * nx + z_low * nx * ny) += x_low_v * y_low_v * z_low_v;
                    //Grid_weight(x_low + y_low * nx + z_high * nx * ny) += x_low_v * y_low_v * z_high_v;
                }
                if (y_high < ny){
                    Grid(x_low + y_high * nx + z_low * nx * ny) += x_low_v * y_high_v * z_low_v * V(row);
                    Grid(x_low + y_high * nx + z_high * nx * ny) += x_low_v * y_high_v * z_high_v * V(row);
                    
                    //Grid_weight(x_low + y_high * nx + z_low * nx * ny) += x_low_v * y_high_v * z_low_v;
                    //Grid_weight(x_low + y_high * nx + z_high * nx * ny) += x_low_v * y_high_v * z_high_v;
                }
            }
            
            if (x_high < nx) {
                if (y_low >= 0){
                    Grid(x_high + y_low * nx + z_low * nx * ny) += x_high_v * y_low_v * z_low_v * V(row);
                    Grid(x_high + y_low * nx + z_high * nx * ny) += x_high_v * y_low_v * z_high_v * V(row);
                    
                    //Grid_weight(x_high + y_low * nx + z_low * nx * ny) += x_high_v * y_low_v * z_low_v;
                    //Grid_weight(x_high + y_low * nx + z_high * nx * ny) += x_high_v * y_low_v * z_high_v;
                }
                if (y_high < ny){
                    Grid(x_high + y_high * nx + z_low * nx * ny) += x_high_v * y_high_v * z_low_v * V(row);
                    Grid(x_high + y_high * nx + z_high * nx * ny) += x_high_v * y_high_v * z_high_v * V(row);
                    
                    //Grid_weight(x_high + y_high * nx + z_low * nx * ny) += x_high_v * y_high_v * z_low_v;
                    //Grid_weight(x_high + y_high * nx + z_high * nx * ny) += x_high_v * y_high_v * z_high_v;
                    
                }
            }
        }
        
    }
    
    //for (int i = 0; i < nx * ny * nz; i ++) {
        //if (Grid_weight(i) != 0) {
            //Grid(i) = Grid(i) / Grid_weight(i);
        //}
    //}
    

}

