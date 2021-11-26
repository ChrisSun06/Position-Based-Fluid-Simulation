#ifndef COEFS_H
#define COEFS_H

struct Coef {
	double MASS; // constant particle mass
	double STIFFNESS; // stiffness coeffcient
	double RHO_IDEAL; // ideal gas rest pressure
	double VISCOSITY; // viscosity coefficient
	double H; // support radius
	double G;
	double SUFRACE_THRESH; // for numerical stability
	double SUFRACE_TENSION; // surface tension
};

#endif