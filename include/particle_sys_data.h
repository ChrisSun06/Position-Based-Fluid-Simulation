#ifndef PARTICLES_DATA_H
#define PARTICLES_DATA_H

#include <Eigen/Dense>

struct Particles {
	Eigen::MatrixXd PredictedPos;
	Eigen::MatrixXd deltaP;
	Eigen::MatrixXd position;
	Eigen::MatrixXd velocity;
	Eigen::MatrixXd acceleration;
	Eigen::MatrixXd force;
	Eigen::VectorXd pressure;
	Eigen::VectorXd density;
	Eigen::VectorXd C;
	Eigen::VectorXd grad_C;
	Eigen::MatrixXd accum_Grad_C;
	Eigen::VectorXd lambda;
};

#endif