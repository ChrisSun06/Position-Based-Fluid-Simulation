#include <vector>
#include <Kernels.h>
#include <PBF_solver.h>
#include <Particles.h>

#define CONST_INV_REST_DENSITY 0.0001
#define RELAXATION 0.01
#define K 0.0003
#define PRESSURE_N 6
#define EPSILON 0.0000325f
#define VISCOSITY_C 0.000045
#define RADIUS 0.07 
#define GRAVITY 9.8
#define MASS 1.0

using namespace std;

void collision_response_2(
	Particles &particles
) {
	Eigen::Vector3d pt;
	Eigen::Vector3d normal;
	double dot_product;
	for (int i = 0; i < particles.q.rows(); i++) {
		Eigen::Vector3d pi = particles.q_new.row(i);
		// wall 1
		pt <<  1.,0.,0.;
		normal << -1., 0., 0.;
		dot_product = (pi - pt).dot(normal);
		if (dot_product < 0.) {
            particles.q_new.row(i) -= 1.1 * dot_product * normal;
			particles.v.row(i) = Eigen::Vector3d::Zero();
        }
		// wall 2
		pt <<  -1.,0.,0.;
		normal <<  1.,0.,0.;
		dot_product = (pi - pt).dot(normal);
		if (dot_product < 0.) {
            particles.q_new.row(i) -= 1.1 * dot_product * normal;
			particles.v.row(i) = Eigen::Vector3d::Zero();
        }
		// wall 3
		pt <<  0.,1.,0.;
		normal << 0.,-1.,0.;
		dot_product = (pi - pt).dot(normal);
		if (dot_product < 0.) {
            particles.q_new.row(i) -= 1.1 * dot_product * normal;
			particles.v.row(i) = Eigen::Vector3d::Zero();
        }
		// wall 4
		pt <<  0.,-1.,0.;
		normal << 0.,1.,0.;
		dot_product = (pi - pt).dot(normal);
		if (dot_product < 0.) {
            particles.q_new.row(i) -= 1.1 * dot_product * normal;
			particles.v.row(i) = Eigen::Vector3d::Zero();
        }
		// wall 5
		pt <<  0.,0.,1.;
		normal << 0.,0.,-1.;
		dot_product = (pi - pt).dot(normal);
		if (dot_product < 0.) {
            particles.q_new.row(i) -= 1.1 * dot_product * normal;
			particles.v.row(i) = Eigen::Vector3d::Zero();
        }
		// wall 6
		pt <<  0.,0.,-1.;
		normal << 0.,0.,1.;
		dot_product = (pi - pt).dot(normal);
		if (dot_product < 0.) {
            particles.q_new.row(i) -= 1.1 * dot_product * normal;
			particles.v.row(i) = Eigen::Vector3d::Zero();
        }
	}
}

void resolve_neighbors(Particles &particles, std::vector< std::vector<int> > &neighbors){
	Eigen::Vector3d maxs = particles.q_new.colwise().maxCoeff();
	Eigen::Vector3d mins = particles.q_new.colwise().minCoeff();
	Eigen::Vector3d ones;
	ones << 1,1,1;
	Eigen::VectorXd step = (maxs - mins) / RADIUS + ones;
	vector<vector<vector<vector<int>>>> neighborhood(int(step(0)), vector<vector<vector<int>>>(int(step(1)), vector<vector<int>>(int(step(2)), vector<int>())));
	std::vector<int> neighbor_particle = {-1,0,1};
	int rows = particles.q.rows();
	Eigen::MatrixXd mins_mat = mins.transpose().replicate(rows, 1);
	Eigen::MatrixXd diff = (particles.q_new - mins_mat) / RADIUS;
	for(int i = 0; i < particles.q.rows(); i++){
		neighborhood[int(diff(i, 0))][int(diff(i, 1))][int(diff(i, 2))].push_back(i);
	}
	for(int i = 0; i < particles.q.rows(); i++){
		for (int offset_x : neighbor_particle) {
			int neighbor_x = offset_x + int(diff(i, 0));
			if (neighbor_x < 0 || neighbor_x > int(step(0))-1){
				continue;
			}
			for (int offset_y : neighbor_particle) {
				int neighbor_y = offset_y + int(diff(i, 1));
				if (neighbor_y < 0 || neighbor_y > int(step(1))-1){
					continue;
				}
				for (int offset_z : neighbor_particle) {
					int neighbor_z = offset_z + int(diff(i, 2));
					if (neighbor_z < 0 || neighbor_z > int(step(2))-1){
						continue;
					}
					for (int b : neighborhood[neighbor_x][neighbor_y][neighbor_z]) {
						if(i != b){
							neighbors[i].push_back(b);
						}
					}
				}
			}
		}
	}
}


void addExtForce(double dt, const Eigen::Vector3d & extF, Particles& particles){
	for (int i = 0; i < particles.q.rows(); i++) {
		particles.v.row(i) += dt*(particles.f.row(i)+extF.transpose())/MASS;
		particles.deltaP.row(i).setZero();
		particles.q_new.row(i) = particles.q.row(i)+particles.v.row(i)*dt;
		particles.f.row(i).setZero();
	}
	collision_response_2(particles);
}


void solveLambda(Particles& particles, std::vector< std::vector<int> > neighbors){
	Eigen::Vector3d p;
	Eigen::Vector3d p_j;
	Eigen::Vector3d g;
	// Solve for C_i
	for (int i = 0; i < particles.q.rows(); i++) {
		double density_i = 0;
		p = particles.q_new.row(i);
		for (int j : neighbors[i])
		{
			Eigen::Vector3d diff;
			p_j = particles.q_new.row(j);
			diff = p-p_j;
			double W = kernel(diff, RADIUS);
			density_i += W;
		}
		particles.density(i) = density_i;
	}

	for (int i = 0; i < particles.q.rows(); i++) {
		double densityConstraint = (particles.density(i) * CONST_INV_REST_DENSITY) - 1.0;
		densityConstraint = densityConstraint > 0 ? densityConstraint : 0.0;
		particles.C(i) = densityConstraint;
	}
	
	// Solve for gradient C_i
	for (int i = 0; i < particles.q.rows(); i++) {
		Eigen::Vector3d grad;
		Eigen::Vector3d diff;
		Eigen::Vector3d grad_k_ci;
		p = particles.q_new.row(i);
		grad.setZero();
		double denom = 0;
		for (int j : neighbors[i])
		{
			p_j = particles.q_new.row(j);
			Eigen::Vector3d grad_spiky;
			grad_spiky.setZero();
			spiky(p-p_j, RADIUS, grad_spiky);
			grad += grad_spiky;
			
		}
		grad *= CONST_INV_REST_DENSITY;
		denom += std::pow(grad.norm(), 2.0); 
		grad.setZero();
		for(int j: neighbors[i]){
			p_j = particles.q_new.row(j);
			spiky(p - p_j, RADIUS, grad);
			grad *= -CONST_INV_REST_DENSITY;
			denom += std::pow(grad.norm(), 2.0);
		}
		double sum_Ci = denom + RELAXATION;

		// Solve for lambda_i
		particles.lambda(i)=-1.0*(particles.C(i)/(sum_Ci));
	}
}

void solvePosition(Particles& particles, std::vector< std::vector<int> > neighbors){
	Eigen::Vector3d p;
	Eigen::Vector3d p_j;
	Eigen::Vector3d delta_p;
	Eigen::Vector3d grad;
	for (int i = 0; i < particles.q.rows(); i++) {
		delta_p.setZero();
		p = particles.q_new.row(i);
		double s_corr=0.0;
		double lambda_i_sum = 0.0;
		double lambda_j_sum = 0.0;
		double s_corr_sum = 0.0;
		for (int j : neighbors[i]) {
			if (i==j) continue;
			p_j = particles.q_new.row(j);
			s_corr = -K * pow(kernel(p-p_j, RADIUS) / Wpoly(0.1*RADIUS, RADIUS), 4.0);
			spiky(p-p_j, RADIUS, grad);
			if(abs(s_corr) < 0.000001){
				s_corr = 0.0;
			}
			double a = particles.lambda(i) + particles.lambda(j) + s_corr;
			delta_p += a * grad;
		}
			
		particles.deltaP.row(i) = CONST_INV_REST_DENSITY * delta_p;
	}
	for (int i = 0; i < particles.q.rows(); i++) {
		particles.q_new.row(i)+=particles.deltaP.row(i);
		particles.deltaP.row(i).setZero();
	}
	collision_response_2(particles);
}


void add_vorticity_ch(Particles& particles, std::vector< std::vector<int> > neighbors){
	Eigen::Vector3d dx;
	Eigen::Vector3d dy;
	Eigen::Vector3d dz;
	double DELTA = 1e-2;
	dx << DELTA, 0, 0;
	dy << 0, DELTA, 0;
	dz << 0, 0, DELTA;
	for (int i = 0; i < particles.q.rows(); i++) {
      	Eigen::Vector3d pi = particles.q_new.row(i);
		Eigen::Vector3d vi = particles.v.row(i);
      	Eigen::Vector3d wi;
		wi.setZero();
      
		Eigen::Vector3d omega_dx;
		Eigen::Vector3d omega_dy;
		Eigen::Vector3d omega_dz;
		omega_dx.setZero();
		omega_dy.setZero();
		omega_dz.setZero();

		for (int j: neighbors[i]) {
			Eigen::Vector3d vj = particles.v.row(j);
			Eigen::Vector3d pj = particles.q_new.row(j);

			Eigen::Vector3d grad_spiky;
			grad_spiky.setZero();
			spiky(pj - pi, RADIUS, grad_spiky);
			
			Eigen::Vector3d vij = vj - vi;
			
			Eigen::Vector3d crossij = grad_spiky.cross(vij);
			
			wi += crossij;

			Eigen::Vector3d pij_dx;
			Eigen::Vector3d pij_dy;
			Eigen::Vector3d pij_dz;

			pij_dx = pj-pi + dx;
			pij_dy = pj-pi + dy;
			pij_dz = pj-pi + dz;

			Eigen::Vector3d grad_x;
			Eigen::Vector3d grad_y;
			Eigen::Vector3d grad_z;
			spiky(pij_dx, RADIUS, grad_x);
			spiky(pij_dy, RADIUS, grad_y);
			spiky(pij_dz, RADIUS, grad_z);

			crossij = grad_x.cross(vij);
			omega_dx += crossij;
			crossij = grad_y.cross(vij);
			omega_dy += crossij;
			crossij = grad_z.cross(vij);
			omega_dz += crossij;
		}
		if (wi.norm() == 0.0) {
			return;
		}
		Eigen::Vector3d eta;
		eta.setZero();
		eta << (omega_dx.norm() - wi.norm()) / DELTA, (omega_dy.norm() - wi.norm()) / DELTA, (omega_dz.norm() - wi.norm()) / DELTA;
		
		eta.normalize();
		
		Eigen::Vector3d f_vorticity;
		f_vorticity = eta.cross(wi);
		f_vorticity *= EPSILON;
		particles.f.row(i) = f_vorticity;
	}
}

void calculate_viscosity(Particles& particles, std::vector< std::vector<int> > neighbors){
	for(int i = 0; i < particles.q.rows(); i++){
		Eigen::Vector3d pi = particles.q_new.row(i);
		Eigen::Vector3d vi = particles.v.row(i);
		Eigen::Vector3d viscosity;
		viscosity.setZero();
		for(int j : neighbors[i]){
			Eigen::Vector3d pj = particles.q_new.row(j);
			Eigen::Vector3d vj = particles.v.row(j);
			Eigen::Vector3d vij = vj - vi;
			viscosity += vij * kernel(pi - pj, RADIUS);
		}
		particles.v.row(i) += VISCOSITY_C*viscosity;
	}
}

void final_update(Particles& particles, double dt, std::vector< std::vector<int> > neighbors){
	for (int i = 0; i < particles.q.rows(); i++) {
		particles.v.row(i)=(particles.q_new.row(i)-particles.q.row(i))/dt;
	}

	calculate_viscosity(particles, neighbors);
	add_vorticity_ch(particles, neighbors);

	for (int i = 0; i < particles.q.rows(); i++) {
		particles.q.row(i)=particles.q_new.row(i);
	}
}

void cleanup(Particles& particles){
	particles.q_new.setZero();
	particles.deltaP.setZero();
	particles.C.setZero();
	particles.grad_C.setZero();
	particles.lambda.setZero();
}

void PBF_update(Particles& particles, double dt) {
	
	std::vector<int> v_row;
	std::vector< std::vector<int> > neighbors(particles.q.rows(), v_row);

	Eigen::Vector3d gravity_force;
	gravity_force  << 0, -GRAVITY * MASS, 0;

	cleanup(particles);

	addExtForce(dt, gravity_force, particles);

	resolve_neighbors(particles, neighbors);

	for(unsigned int iter=0;iter<3;iter++){
		solveLambda(particles, neighbors);
		solvePosition(particles, neighbors);
	}
	final_update(particles, dt, neighbors);
}
