#include "initialize_particles.h"

void initialize_particles(
    Eigen::MatrixXd & P,
    const Eigen::RowVector3d Starting_corner,
    const double step,
    const int nx,
    const int ny,
    const int nz)
{
    P.resize(nx * ny * nz, 3);
    Eigen::RowVector3d corner = Starting_corner;
    int counter = 0;
    
    for (int i = 0; i < nx; i ++){
        corner(1) = Starting_corner(1);
        for (int j = 0; j < ny; j ++){
            corner(2) = Starting_corner(2);
            for (int k = 0; k < nz; k ++){
                P.row(counter) = corner;
                counter += 1;
                corner(2) += step;
            }
            corner(1) += step;
        }
        corner(0) += step;
    }
}
