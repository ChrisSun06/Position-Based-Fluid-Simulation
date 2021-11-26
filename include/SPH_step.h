#ifndef SPH_STEP_H
#define SPH_STEP_H

#include <Eigen/Dense>
#include <particle_sys_data.h>
#include <coefs.h>
#include <vector>

// void SPH_step(Particles& particles, double dt, Coef& coef);
void PBF_update(Particles& particles, double dt, Coef& coef);
typedef std::vector<std::vector<std::vector<int>>> Grid;

#endif
