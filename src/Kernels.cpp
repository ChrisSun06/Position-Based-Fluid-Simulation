#include <Kernels.h>
#include <iostream>

double kernel(Eigen::Ref<const Eigen::Vector3d> r, double H){
	if(r.norm()>H){
		return 0.00000;
	}
	return 315. / (64.*M_PI*pow(H, 9.0)) * pow(H*H - r.dot(r), 3.0);
}

void spiky(
	Eigen::Ref<const Eigen::Vector3d> r, double H, Eigen::Vector3d& g
) {
	if(r.norm()>H){
		Eigen::Vector3d ret;
		ret << 0.000000, 0.000000, 0.000000;
		g = ret;
		return;
	}
	g = -45. / (M_PI*pow(H, 6.0)) * pow(H - r.norm(), 2.0) * r/(r.norm()+1e-9);
}

double Wpoly(double r, double H) {
	if(r>H){
		return 0.00000;
	}
	return 315. / (64. * M_PI * pow(H, 9.0)) * pow(H*H - r*r, 3.0);
}