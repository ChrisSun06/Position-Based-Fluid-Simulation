#include <SPH_step.h>
#include <vector>
#include <weight_funcs.h>
#include <coefs.h>
#include <iostream>

#define CONST_INV_REST_DENSITY 0.0001
#define RELAXATION 0.01
#define PRESSURE_K 0.000001
#define PRESSURE_N 6
#define EPSILON 0.00005f
#define VISCOSITY_C 0.00001

Eigen::MatrixXd N(6, 3);
Eigen::MatrixXd P(6, 3);

double
	density,
	d2Cs, // Color surface Laplacian
	r; // positions distance
Eigen::Vector3d d,
	f_pressure, // pressure 
	f_visco, // viscosity 
	f_surface, // surface tension 
	f_g, // gravity 
	grad_press, // W_pressure gradient
	grad_poly, // W_poly gradient
	dCs; // Color surface normal 

double min_x;
double min_y;
double min_z;
double max_x;
double max_y;
double max_z;

int cell_x;
int cell_y;
int cell_z;

std::vector<std::vector<std::vector<std::vector<int> > > > grid(cell_x, std::vector<std::vector<std::vector<int> > >(cell_y, std::vector<std::vector<int> >(cell_z, std::vector<int>())));
int grid_x, grid_y, grid_z;

void collision_response_2(
	Particles &particles,
	Eigen::Ref<const Eigen::MatrixXd> P,
	Eigen::Ref<const Eigen::MatrixXd> N
) {
	for (int i = 0; i < particles.position.rows(); i++) {
		for (int w = 0; w < P.rows(); w++) {
			Eigen::Vector3d pi = particles.PredictedPos.row(i);
			Eigen::Vector3d vi = particles.velocity.row(i);
			Eigen::Vector3d p_wall = P.row(w);
			Eigen::Vector3d N_wall = N.row(w);
			double d = (p_wall - pi).dot(N_wall);
			Eigen::Vector3d proj;

			if (d > 0.) {
				// double nom = (N_wall(0) * p_wall(0)) - (N_wall(0) * pi(0)) + (N_wall(1) * p_wall(1)) - (N_wall(1) * pi(1)) + (N_wall(2) * p_wall(2)) - (N_wall(2) * pi(2));
				// double denom = pow(N_wall(0), 2) + pow(N_wall(1), 2) + pow(N_wall(2), 2);
				// double t = nom / denom;
				// proj = pi + t * 2 * N_wall;
				// push particle out of wall
				particles.PredictedPos.row(i) += 1.1 * d * N_wall;
				particles.velocity.row(i) = Eigen::Vector3d::Zero();
				// particles.PredictedPos.row(i) = proj;
				// reflect the velocity component
				// particles.velocity.row(i) -=  1.2 * vi.dot(N_wall) * N_wall;
				// particles.velocity.row(i) -=  1.0*vi.dot(N_wall) * N_wall;
				// particles.PredictedPos.row(i) += particles.velocity.row(i) * 0.02;
			}
		}
	}

}


void addExtForce(double dt, const Eigen::Vector3d & extF, Particles& particles, Coef& coef){
	N << 0., -1., 0.,
	0., 1., 0.,
	1., 0., 0.,
	0., 0., 1.,
	-1., 0., 0.,
	0., 0., -1.;
	P <<
	0., 3., 0.,
	0., -1., 0.,
	-1., 0., 0.,
	0, 0., -1.,
	1., 0., 0.,
	0., 0., 1.;
	for (int i = 0; i < particles.position.rows(); i++) {
		particles.velocity.row(i) += dt*(particles.force.row(i)+extF.transpose())/coef.MASS;
		particles.deltaP.row(i).setZero();
		particles.PredictedPos.row(i) = particles.position.row(i)+particles.velocity.row(i)*dt;
		particles.force.row(i).setZero();
	}
	collision_response_2(particles, P, N);
}

void findNeighbors(std::vector< std::vector<int> > &neighbors){
	for (grid_x = 0; grid_x < cell_x; grid_x++) {
		for (grid_y = 0; grid_y < cell_y; grid_y++) {
			for (grid_z = 0; grid_z < cell_z; grid_z++) {
				assert(grid_x < grid.size());
				assert(grid_y < grid[grid_x].size());
				assert(grid_z < grid[grid_x ][grid_y ].size());
				for (int i : grid[grid_x][grid_y][grid_z]) {
					// neighbors.push_back(std::vector<int>());
					for (int count_x = -1; count_x < 2; count_x++) {
						if (grid_x + count_x < 0 || grid_x + count_x >= cell_x) continue;
						for (int count_y = -1; count_y < 2; count_y++) {
							if (grid_y + count_y < 0 || grid_y + count_y >= cell_y) continue;
							for (int count_z = -1; count_z < 2; count_z++) {
								if (grid_z + count_z < 0 || grid_z + count_z >= cell_z) continue;
								assert(grid_x + count_x >= 0 && grid_x + count_x < grid.size());
								assert(grid_y + count_y >= 0 && grid_y + count_y < grid[grid_x + count_x ].size());
								assert(grid_z + count_z >= 0 && grid_z + count_z < grid[grid_x + count_x][grid_y + count_y].size());
								for (int j : grid[grid_x + count_x][grid_y + count_y][grid_z + count_z]) {
									if(i != j){
										neighbors[i].push_back(j);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}



void solveLambda(Particles& particles, Coef& coef, std::vector< std::vector<int> > neighbors){
	Eigen::Vector3d p;
	Eigen::Vector3d p_j;
	Eigen::Vector3d g;
	/////////////////Solve C_i///////////////////////////////////////////////
	for (int i = 0; i < particles.position.rows(); i++) {
		// p = particles.PredictedPos.row(i);
		// double ro=0.0;
		// for(int j : neighbors[i]){
		// 		if(i == j) continue;
		// 		p_j = particles.PredictedPos.row(j);
		// 		double rr = sqrt((p-p_j).dot(p-p_j));
		// 		// double kv=kernel(p-p_j,coef.H);
		// 		double kv = Wpoly(rr, coef.H);
		// 		particles.density(i) += kv * coef.MASS;
		// 		// particles.density(j) += kv * coef.MASS;
		// }
		// Eigen::Vector3d tmp;
		// tmp.setZero();
		// particles.density(i) += kernel(tmp,coef.H) * coef.MASS;
		// Particles[i]->density+=kernel(vec3(),core_radius)/Particles[i]->inv_mass;
		double density_i = 0;
		p = particles.PredictedPos.row(i);
		for (int j : neighbors[i])
		{
			//if(i == j) continue;
			Eigen::Vector3d diff;
			p_j = particles.PredictedPos.row(j);
			diff = p-p_j;
			double W = kernel(diff, coef.H);
			// double mag = sqrt(diff.transpose()*diff);
			// //calculate poly6 kernal, assume all particles have unit mass
			// double W = 315. / (64 * 3.14*pow(coef.H, 9));
			// if (mag >= 0 && mag <= coef.H)
			// {
			// 	W *= pow((pow(coef.H, 2) - pow(mag, 2)), 3);
			// }
			// else
			// {
			// 	W *= 0;
			// }
			density_i += W;
		}
		particles.density(i) = density_i;
	}

	for (int i = 0; i < particles.position.rows(); i++) {
		double densityConstraint = (particles.density(i) * CONST_INV_REST_DENSITY) - 1.0;
		// printf("lambda_C%f\n", densityConstraint);
		densityConstraint = densityConstraint > 0 ? densityConstraint : 0.0;
		// particles.C(i)=particles.density(i)*CONST_INV_REST_DENSITY-1.0;
		// particles.C(i) = particles.C(i) > 0 ? particles.C(i) : 0;
		// particles.density(i)=0.0;
		particles.C(i) = densityConstraint;
	}
	
	/////////////////////////////////////////////////////////////////////////

	////////////////////////Solve C_i_gradient///////////////////////////////
	
	// for (int i = 0; i < particles.position.rows(); i++) {
	// 	p = particles.PredictedPos.row(i);
	// 	double C_i_gradient, sum_gradients=0.0;
	// 	for(int j : neighbors[i]){
				
	// 		if(i == j) continue;
	// 		p_j = particles.PredictedPos.row(j);
	// 		spiky(p-p_j,coef.H, g);
	// 		C_i_gradient=std::pow(g.norm()*CONST_INV_REST_DENSITY,2);
	// 		//sum_gradients+=C_i_gradient;
	// 		particles.grad_C(i)+=C_i_gradient;
	// 		// particles.grad_C(j)+=C_i_gradient;
	// 	}

	// 	for(int j : neighbors[i]){
	// 		if(i == j) continue;
	// 		p_j = particles.PredictedPos.row(j);
	// 		Eigen::Vector3d agc;
	// 		spiky(p-p_j,coef.H, agc);
	// 		agc = agc * CONST_INV_REST_DENSITY;
	// 		particles.accum_Grad_C.row(i)+=agc;
	// 		particles.accum_Grad_C.row(j)-=agc;
	// 	}
	// 	C_i_gradient=particles.accum_Grad_C.row(i).norm()*CONST_INV_REST_DENSITY;
	// 	particles.grad_C(i)+=C_i_gradient*C_i_gradient;

	for (int i = 0; i < particles.position.rows(); i++) {
		Eigen::Vector3d grad;
		Eigen::Vector3d diff;
		Eigen::Vector3d grad_k_ci;
		p = particles.PredictedPos.row(i);
		grad.setZero();
		double denom = 0;
		for (int j : neighbors[i])
		{
			// if(i == j) continue;
			p_j = particles.PredictedPos.row(j);
			// diff = p-p_j;
			// double mag = sqrt(diff.transpose()*diff);
			// double W = -45.0 / (3.14*pow(coef.H, 6.0));
			// if (mag >= 0 && mag <= coef.H)
			// {
			// 	// P3D normalized = P3D(diff[0] / mag, diff[1] / mag, diff[2] / mag);
			// 	W *= std::pow(coef.H - mag, 2.0);
			// 	Eigen::Vector3d tmp = 1000*(diff/(diff.norm() + 0.001))*W;
			// 	grad -= tmp;
			// 	double magGrad = sqrt(tmp.transpose()*tmp);
			// 	denom += std::pow(magGrad, 2.0);
			// }
			Eigen::Vector3d grad_spiky;
			grad_spiky.setZero();
			spiky(p-p_j, coef.H, grad_spiky);
			grad += grad_spiky;
			
		}
		grad *= CONST_INV_REST_DENSITY;
		denom += std::pow(grad.norm(), 2.0); 
		grad.setZero();
		for(int j: neighbors[i]){
			p_j = particles.PredictedPos.row(j);
			spiky(p - p_j, coef.H, grad);
			grad *= -CONST_INV_REST_DENSITY;
			denom += std::pow(grad.norm(), 2.0);
		}
		// else {
		// 		Eigen::Vector3d grad_spiky;
		// 		spiky(p-p_j, coef.H, grad);
		// 		grad_spiky = -1.0 * CONST_INV_REST_DENSITY * grad;
		// 		denom += std::pow(grad_spiky.norm(), 2.0);
		// 	}
		// double mG = sqrt(grad.transpose()*grad);
		// double sum_Ci = denom + pow(mG, 2) + RELAXATION;
		double sum_Ci = denom + RELAXATION;
		// printf("lambda_C%f\n", particles.C(i));
		// printf("lambda_Ci%f\n", sum_Ci);
		// double sum_Ci=particles.grad_C(i)+RELAXATION;
		particles.lambda(i)=-1.0*(particles.C(i)/(sum_Ci));
		// if(particles.lambda(i) != 0.0){
		// 	std::cout << "Ci: " << particles.C(i) << " sum_ci: " << sum_Ci << "\n";
		// }
		// particles.grad_C(i)=0.0;
		// particles.accum_Grad_C.row(i).setZero();
	}
}

void solvePosition(Particles& particles, Coef& coef, std::vector< std::vector<int> > neighbors){
	Eigen::Vector3d p;
	Eigen::Vector3d p_j;
	Eigen::Vector3d delta_p;
	Eigen::Vector3d grad;
	for (int i = 0; i < particles.position.rows(); i++) {
		delta_p.setZero();
		p = particles.PredictedPos.row(i);
		// Eigen::Vector3d tmp;
		// tmp << 1.0, 1.0, 1.0;
		// Eigen::Vector3d d_q= (0.1*coef.H)*tmp;
		double s_corr=0.0;
		double lambda_i_sum = 0.0;
		double lambda_j_sum = 0.0;
		double s_corr_sum = 0.0;
		for (int j : neighbors[i]) {
			if (i==j) continue;
			p_j = particles.PredictedPos.row(j);
			s_corr = -PRESSURE_K * pow(kernel(p-p_j, coef.H) / Wpoly(0.1*coef.H, coef.H), 4.0);
			spiky(p-p_j, coef.H, grad);
			if(abs(s_corr) < 0.000001){
				s_corr = 0.0;
			}
			double a = particles.lambda(i) + particles.lambda(j) + s_corr;
			// if(particles.lambda(i) != 0.0 && particles.lambda(j) != 0.0 && s_corr != 0.0){
			// 	printf("i: %f, j: %f, S_corr: %f\n", particles.lambda(i), particles.lambda(j), s_corr);
			// 	std::cout << "a: " << a << " Grad: " << grad << "\n";
			// }
			// printf("lambda%f\n", particles.lambda(i));
			delta_p += a * grad;
		}
			
		particles.deltaP.row(i) = CONST_INV_REST_DENSITY * delta_p;
	}
	for (int i = 0; i < particles.position.rows(); i++) {
		particles.PredictedPos.row(i)+=particles.deltaP.row(i);
		particles.deltaP.row(i).setZero();
	}
	N << 0., -1., 0.,
	0., 1., 0.,
	1., 0., 0.,
	0., 0., 1.,
	-1., 0., 0.,
	0., 0., -1.;
	P <<
	0., 3., 0.,
	0., -1., 0.,
	-1., 0., 0.,
	0, 0., -1.,
	1., 0., 0.,
	0., 0., 1.;
	collision_response_2(particles, P, N);
}

void calculate_vorticity(Particles& particles, Coef& coef, std::vector< std::vector<int> > neighbors){
	for(int i = 0; i < particles.position.rows(); i++){
		
		Eigen::Vector3d pi = particles.PredictedPos.row(i);
		Eigen::Vector3d vi = particles.velocity.row(i);
		Eigen::Vector3d wi;
		wi.setZero();
		for(int j: neighbors[i]){
			Eigen::Vector3d vj = particles.velocity.row(j);
			Eigen::Vector3d pj = particles.PredictedPos.row(j);
			Eigen::Vector3d grad_spiky;
			grad_spiky.setZero();
			spiky(pi - pj, coef.H, grad_spiky);
			wi += (vj - vi).cross(grad_spiky);
		}
		Eigen::Vector3d eta;
		eta.setZero();
		// wi.normalized();
		for(int j: neighbors[i]){
			Eigen::Vector3d pj = particles.PredictedPos.row(j);
			Eigen::Vector3d grad_spiky;
			grad_spiky.setZero();
			spiky(pi - pj, coef.H, grad_spiky);
			eta += grad_spiky*wi.norm();
		}
		if(eta.norm() == 0.0){
			return;
		}
		else{
			Eigen::Vector3d normalized = eta / eta.norm();
			Eigen::Vector3d result = normalized.cross(wi);
			particles.force.row(i) += EPSILON * result;
		}
 		
	}
}

void calculate_viscosity(Particles& particles, Coef& coef, std::vector< std::vector<int> > neighbors){
	for(int i = 0; i < particles.position.rows(); i++){
		Eigen::Vector3d pi = particles.PredictedPos.row(i);
		Eigen::Vector3d vi = particles.velocity.row(i);
		Eigen::Vector3d viscosity;
		viscosity.setZero();
		for(int j : neighbors[i]){
			Eigen::Vector3d pj = particles.PredictedPos.row(j);
			Eigen::Vector3d vj = particles.velocity.row(j);
			Eigen::Vector3d vij = vj - vi;
			viscosity += vij * kernel(pi - pj, coef.H);
		}
		particles.velocity.row(i) += VISCOSITY_C*viscosity;
	}
}
void add_vorticity_ch(Particles& particles, Coef& coef, std::vector< std::vector<int> > neighbors){
	Eigen::Vector3d dx;
	Eigen::Vector3d dy;
	Eigen::Vector3d dz;
	double DELTA = 0.01;
	dx << DELTA, 0, 0;
	dy << 0, DELTA, 0;
	dz << 0, 0, DELTA;
	for (int i = 0; i < particles.position.rows(); i++) {
      	Eigen::Vector3d pi = particles.PredictedPos.row(i);
		Eigen::Vector3d vi = particles.velocity.row(i);
      	Eigen::Vector3d wi;
		wi.setZero();
      
		Eigen::Vector3d omegaiPurturbedX;
		Eigen::Vector3d omegaiPurturbedY;
		Eigen::Vector3d omegaiPurturbedZ;
		omegaiPurturbedX.setZero();
		omegaiPurturbedY.setZero();
		omegaiPurturbedZ.setZero();

		for (int j: neighbors[i]) {
			Eigen::Vector3d vj = particles.velocity.row(j);
			Eigen::Vector3d pj = particles.PredictedPos.row(j);
			// Eigen::Vector3d pij = new Vector3d(pj.x);
			// pij.sub(pi.x);
			// Eigen::Vector3d dWdpj = Kernels.poly6dp(pij);
			Eigen::Vector3d grad_spiky;
			grad_spiky.setZero();
			spiky(pj - pi, coef.H, grad_spiky);
			
			// Vector3d vij = new Vector3d(pj.v);
			// vij.sub(pi.v);
			Eigen::Vector3d vij = vj - vi;
			
			// Vector3d crossij = new Vector3d();  
			// crossij.cross(dWdpj, vij);
			Eigen::Vector3d crossij = grad_spiky.cross(vij);
			
			// omegai.add(crossij);
			wi += crossij;
			
			// Vector3d pijPurturbedX = new Vector3d(pij);
			// Vector3d pijPurturbedY = new Vector3d(pij);
			// Vector3d pijPurturbedZ = new Vector3d(pij);
			Eigen::Vector3d pijPurturbedX;
			Eigen::Vector3d pijPurturbedY;
			Eigen::Vector3d pijPurturbedZ;
			// pijPurturbedX.add(dx);
			// pijPurturbedY.add(dy);
			// pijPurturbedZ.add(dz);
			pijPurturbedX = pj-pi + dx;
			pijPurturbedY = pj-pi + dy;
			pijPurturbedZ = pj-pi + dz;

			Eigen::Vector3d grad_x;
			Eigen::Vector3d grad_y;
			Eigen::Vector3d grad_z;
			spiky(pijPurturbedX, coef.H, grad_x);
			spiky(pijPurturbedY, coef.H, grad_y);
			spiky(pijPurturbedZ, coef.H, grad_z);

			// crossij = crossij.cross(Kernels.poly6dp(pijPurturbedX), vij);
			// omegaiPurturbedX.add(crossij);
			// crossij.cross(Kernels.poly6dp(pijPurturbedY), vij);
			// omegaiPurturbedY.add(crossij);
			// crossij.cross(Kernels.poly6dp(pijPurturbedZ), vij);
			// omegaiPurturbedZ.add(crossij);
			crossij = grad_x.cross(vij);
			omegaiPurturbedX += crossij;
			crossij = grad_y.cross(vij);
			omegaiPurturbedY += crossij;
			crossij = grad_z.cross(vij);
			omegaiPurturbedZ += crossij;
		}
		if (wi(0) == 0.0 && wi(1)  == 0.0 && wi(2)  == 0.0) {
			return;
		}
		Eigen::Vector3d eta;
		eta.setZero();
		eta << (omegaiPurturbedX.norm() - wi.norm()) / DELTA, (omegaiPurturbedY.norm() - wi.norm()) / DELTA, (omegaiPurturbedZ.norm() - wi.norm()) / DELTA;

		// Vector3d eta = new Vector3d(
		// 	(omegaiPurturbedX.length() - omegai.length()) / DELTA,
		// 	(omegaiPurturbedY.length() - omegai.length()) / DELTA,
		// 	(omegaiPurturbedZ.length() - omegai.length()) / DELTA);
		
		eta.normalize();
		
		Eigen::Vector3d fVorticity;
		// fVorticity.cross(eta, omegai);
		fVorticity = eta.cross(wi);
		// fVorticity.scale(vorticityEps);
		fVorticity *= EPSILON;
		particles.force.row(i) = fVorticity;
	}
}

void add_vorticity(Particles& particles, Coef& coef, std::vector< std::vector<int> > neighbors){
	for(int i = 0; i < particles.position.rows(); i++){
		
		Eigen::Vector3d pi = particles.PredictedPos.row(i);
		Eigen::Vector3d vi = particles.velocity.row(i);
		Eigen::Vector3d wi;
		wi.setZero();
		for(int j: neighbors[i]){
			Eigen::Vector3d vj = particles.velocity.row(j);
			Eigen::Vector3d pj = particles.PredictedPos.row(j);
			Eigen::Vector3d grad_spiky;
			grad_spiky.setZero();
			spiky(pi - pj, coef.H, grad_spiky);
			wi += (vj - vi).cross(grad_spiky);
		}
		Eigen::Vector3d eta;
		eta.setZero();
		// wi.normalized();
		for(int j: neighbors[i]){
			Eigen::Vector3d pj = particles.PredictedPos.row(j);
			Eigen::Vector3d p_mid = (pi+pj) / 2.0;
			eta += p_mid - pi;
		}
		if(eta.norm() == 0.0){
			return;
		}
		else{
			eta.normalize();
			particles.force.row(i) += EPSILON * (eta.cross(wi));
		}
 		
	}
}

void final_update(Particles& particles, double dt, Coef& coef, std::vector< std::vector<int> > neighbors){
	for (int i = 0; i < particles.position.rows(); i++) {
		// if(Particles[i]->isOnSurface)
		// 	continue;
		particles.velocity.row(i)=(particles.PredictedPos.row(i)-particles.position.row(i))/dt;
	}
	// TODO: ADD vorticit
	calculate_viscosity(particles, coef, neighbors);
	add_vorticity_ch(particles, coef, neighbors);
	// TODO: ADD viscosi
	for (int i = 0; i < particles.position.rows(); i++) {
		// if(Particles[i]->isOnSurface)
		// 	continue;
		particles.position.row(i)=particles.PredictedPos.row(i);
	}

	//apply_vorticity(dt);
	//apply_viscosity();
	// for(unsigned int i=0;i<ParticleCont;i++){
	// 	if(!dropToGround&&Particles[i]->isOnSurface)
	// 		continue;
	// 	Particles[i]->position=Particles[i]->PredictedPos;
	// }
}

void cleanup(Particles& particles){
	particles.PredictedPos.setZero();
	particles.deltaP.setZero();
	particles.C.setZero();
	particles.grad_C.setZero();
	particles.accum_Grad_C.setZero();
	particles.lambda.setZero();
	// particles.force.setZero();
}

void PBF_update(Particles& particles, double dt, Coef& coef) {
	
	std::vector<int> v_row;
	std::vector< std::vector<int> > neighbors(particles.position.rows(), v_row);

	f_g << 0, coef.G* 1.0, 0;

	cleanup(particles);

	addExtForce(dt, f_g, particles, coef);

	min_x = particles.PredictedPos.col(0).minCoeff();
	min_y = particles.PredictedPos.col(1).minCoeff();
	min_z = particles.PredictedPos.col(2).minCoeff();


	max_x = particles.PredictedPos.col(0).maxCoeff();
	max_y = particles.PredictedPos.col(1).maxCoeff();
	max_z = particles.PredictedPos.col(2).maxCoeff();

	cell_x = (max_x - min_x) / coef.H + 1;
	cell_y = (max_y - min_y) / coef.H + 1;
	cell_z = (max_z - min_z) / coef.H + 1;

	grid = std::vector<std::vector<std::vector<std::vector<int> > > >(cell_x, std::vector<std::vector<std::vector<int> > >(cell_y, std::vector<std::vector<int> >(cell_z, std::vector<int>())));


	for (size_t i = 0; i < particles.PredictedPos.rows(); i++) {
		Eigen::Vector3d pos = particles.PredictedPos.row(i);
		grid_x = (pos(0) - min_x) / coef.H;
		grid_y = (pos(1) - min_y) / coef.H;

		grid_z = (pos(2) - min_z) / coef.H;
		assert(grid_x < grid.size());
		assert(grid_y < grid[grid_x].size());
		assert(grid_z < grid[grid_x][grid_y].size());
		grid[grid_x][grid_y][grid_z].push_back(i);
	}

	findNeighbors(neighbors);

	for(unsigned int iter=0;iter<3;iter++){
		solveLambda(particles, coef, neighbors);
		solvePosition(particles, coef, neighbors);
	}
	final_update(particles, dt, coef, neighbors);
}
