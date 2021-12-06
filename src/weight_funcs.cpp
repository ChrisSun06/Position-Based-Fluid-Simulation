#include <weight_funcs.h>
#include <iostream>

double kernel(Eigen::Ref<const Eigen::Vector3d> r, double H){
	if(r.norm()>H){
		return 0.00000;
	}
	// printf("asd%f\n", r.norm());
	return 315. / (64. * M_PI * pow(H, 9)) * pow(H * H - r.dot(r), 3);
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
	g = -45. / (M_PI * pow(H, 6)) * pow(H - r.norm(), 2) * r / (r.norm() +0.0001);
}

void dWpress(
	Eigen::Ref<const Eigen::Vector3d> d,
	double r, double H,
	Eigen::Vector3d& g
) {
	g = -45. / (M_PI * pow(H, 6)) * pow(H - r, 2) / r * d;
}


double Wpoly(double r, double H) {
	if(r>H){
		return 0.00000;
	}
	// printf("help\n");
	return 315. / (64. * M_PI * pow(H, 9)) * pow(H * H - r * r, 3);
}

void dWpoly(
	Eigen::Ref<const Eigen::Vector3d> d,
	double r, double H,
	Eigen::Vector3d& g
) {
	g = -945. / (32. * M_PI * pow(H, 9)) * pow(H * H - r * r, 2) * d;
}

double d2Wpoly(double r, double H) {
	return -945. / (32. * M_PI * pow(H, 9)) * (H * H - r * r) * (3. * H * H - 7. * r * r);
}

double d2Wvisco(double r, double H) {
	return 45. / (M_PI * pow(H, 6)) * (H - r);
}