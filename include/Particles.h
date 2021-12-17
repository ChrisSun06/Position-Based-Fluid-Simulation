#ifndef Particle
#define Particle

#include <Eigen/Dense>

struct Particles {
	Eigen::MatrixXd q_new;
	Eigen::MatrixXd deltaP;
	Eigen::MatrixXd q;
	Eigen::MatrixXd v;
	Eigen::MatrixXd f;
	Eigen::VectorXd density;
	Eigen::VectorXd C;
	Eigen::VectorXd grad_C;
	Eigen::VectorXd lambda;
};

#endif