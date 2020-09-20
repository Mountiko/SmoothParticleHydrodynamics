
#include <iostream>
#include <assert.h>
#include <sstream>
#include <fstream>
#include <iterator>
#include <cmath>
#include <ctime>
#include <omp.h>
#include "../includes/SPH_2D.h"

using namespace std;

SPH_main* SPH_particle::main_data;

double SPH_particle::cal_mass(double dx, double rho)
{
	mass = pow(dx, 2.0) * rho;
	return mass;
}

double SPH_particle::cal_pressure(double rho, double gamma, double c0)
{
	P = (1000.0 * pow(c0, 2.0) / gamma) * (pow(rho / 1000, gamma) - 1.0);
	return P;
}
void SPH_particle::navier(SPH_particle* neigh, double dx, double gamma, double c0, double mu)
{
	double e_ij[2];
	double v_ij[2];
	double first = 0; //first term in navier-stoke acceleration equation
	double second = 0; //second term

	//cout << neigh.size() << endl;

	//calculate first term without e_ij
	for (int j = 0; j < 2; j++)
	{
		e_ij[j] = (this->x[j] - neigh->x[j]);
		v_ij[j] = this->v[j] - neigh->v[j];
	}

	double dist = sqrt(e_ij[0] * e_ij[0] + e_ij[1] * e_ij[1]);

	if (dist > dx/100 && neigh->mass > rho * dx * dx / 100)
	{
		double dW_dr = main_data->dW_dr(dist);

		first = neigh->mass * (this->P / pow(this->rho, 2.0) + neigh->P / pow(neigh->rho, 2.0)) * dW_dr;
		//calculate second term without v_ij/|r_ij|
		second = neigh->mass * (1 / pow(this->rho, 2.0) + 1 / pow(neigh->rho, 2.0)) * dW_dr / dist;
		//compute e_ij = (x_i-x_j)/|x_i-x_j| and v_ij

		

		//|r_ij| = |x_i-x_j| = dist

		for (int j = 0; j < 2; j++)
		{
			this->a[j] += -first * e_ij[j] / dist + mu * second * v_ij[j];
		}

		D += neigh->mass * dW_dr * (v_ij[0] * e_ij[0] + v_ij[1] * e_ij[1]) / dist;
	}

	// this->a[0] = -first_sum[0] + mu * second_sum[0];
	// this->a[1] = -first_sum[1] + mu * second_sum[1] - 9.81;

}

/*ouble SPH_particle::navier_rho(vector<SPH_particle*> neigh, double dx)

{
	double e_ij[2];
	double v_ij[2];
	double coeffe; //the part without v_ij . e_ij
	double sum = 0;
	double dot = 0; //dot product of v_ij and e_ij
	for (int i = 0; i < neigh.size(); i++)
	{
		coeffe = neigh[i]->cal_mass(dx, neigh[i]->rho) * dW_dr();
		//calculate v_ij and e_ij
		for (int j = 0; j < 2; j++)
		{
			e_ij[j] = this->x[j] - neigh[i]->x[j];
			v_ij[j] = this->v[j] - neigh[i]->v[j];
		}

		double dist = sqrt(e_ij[0] * e_ij[0] + e_ij[1] * e_ij[1]);
		for (int j = 0; j < 2; j++)
			e_ij[j] = e_ij[j] / dist;

		//calculate dot_product
		for (int i = 0; i < 2; i++)
			dot += v_ij[i] * e_ij[i];

		//add to the sum
		sum += coeffe * dot;
	}
	return sum;

}*/


SPH_main::SPH_main(SPH_main* domain)
{
	c0 = domain->c0;
    gamma = domain->gamma;

    mu = domain->mu; // viscosity

    sim_time = domain->sim_time; // time of simulation



    h_fac = domain->h_fac;
    dx = domain->dx;                                //particle initial spacing

    min_x[0] = domain->min_x[0]; 
    min_x[1] = domain->min_x[1]; 
    max_x[0] = domain->max_x[0]; 
    max_x[1] = domain->max_x[1]; 
   //dimensions of simulation region



    max_list[0] = domain->max_list[0]; 
    max_list[1] = domain->max_list[1];

    vector<SPH_particle> particle_list;				//list of all the particles

    vector<vector<vector<SPH_particle*> > > search_grid;
}

SPH_main::SPH_main()
{
	SPH_particle::main_data = this;

	// open the data file
	ifstream bound("./BC");
	ifstream init("./IC");
	ifstream system("./system");

	pair<double, double> init_velocity;
	double rho_in;

	// check if data file is open
	if (bound.is_open() && init.is_open() && system.is_open())
	{
		// set string to get loop through lines of data-file
		string line;
		// Boundary conditions
		bool boundary = false;
		bool is_dx = false;
		bool is_mu = false;
		bool speed_sound = false;
		// Initial conditions
		bool water = false;
		bool velocity = false;
		bool is_rho = false;
		// System variabes
		bool is_h_fac = false;
		bool is_time = false;
		bool is_gamma = false;
		bool is_timestep = false;
		bool is_time_meth = false;
		bool is_adapt_time = false;
		// add a IC and BC count to keep track of what is being captured ------------------ IMPORTANT		
		// loop through the lines of the system-file
		while (getline(system, line))
		{
			if (line.find("#", 0) != string::npos)
			{
				is_h_fac = false;
				is_time = false;
				is_gamma = false;
				is_timestep = false;
				is_time_meth = false;
				is_adapt_time = false;
			}


			if (is_h_fac) h_fac = atof(line.c_str());
			else if (is_time) this->sim_time = atof(line.c_str());
			else if (is_gamma) this->gamma = atof(line.c_str());
			else if (is_timestep) this->dt = atof(line.c_str());
			else if (is_time_meth) this->time_step_meth = line.c_str();
			else if (is_adapt_time)
				istringstream(line.c_str()) >> boolalpha >> this->adaptive_timestepping;

			if (line.find("h_fac:", 0) != string::npos)
				is_h_fac = true;
			else if (line.find("end time[s]:", 0) != string::npos)
				is_time = true;
			else if (line.find("gamma[-]:", 0) != string::npos)
				is_gamma = true;
			else if (line.find("timestep[s]:", 0) != string::npos)
				is_timestep = true;
			else if (line.find("time stepping method:", 0) != string::npos)
				is_time_meth = true;
			else if (line.find("adaptive time stepping:", 0) != string::npos)
				is_adapt_time = true;
		}
		
		// loop through the lines of the BC-file
		while (getline(bound, line))
		{
			// enter this loop once flag has been flipped
			// inside the "if"-loop the relevant data will be extracted

			if (line.find("#", 0) != string::npos)
			{
				boundary = false;
				is_dx = false;
				is_mu = false;
				speed_sound = false;
			}

			if (boundary)
			{
				// tokenize the line to split it into single string
				// of data
				istringstream iss(line);
				vector<string> tokens;
				copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter(tokens));

				pair<double, double> tmp_coord;

				tmp_coord.first = atof(tokens[0].c_str());
				tmp_coord.second = atof(tokens[1].c_str());

				this->boundaryCoords.push_back(tmp_coord);
			}
			else if (is_dx) dx = atof(line.c_str());
			else if (is_mu) this->mu = atof(line.c_str());
			else if (speed_sound) this->c0 = atof(line.c_str());

			if (line.find("boundary[m]:", 0) != string::npos)
				boundary = true;
			else if (line.find("dx[m]:", 0) != string::npos)
				is_dx = true;
			else if (line.find("viscosity[Ns/m2]:", 0) != string::npos)
				is_mu = true;
			else if (line.find("speed of sound[m/s]:", 0) != string::npos)
				speed_sound = true;
		}
		
		// loop through the lines of the IC-file
		while (getline(init, line))
		{
			// enter this loop once flag has been flipped
			// inside the "if"-loop the relevant data will be extracted
			if (line.find("#", 0) != string::npos)
			{
				water = false;
				velocity = false;
				is_rho = false;
			}

			if (water)
			{
				// tokenize the line to split it into single string
				// of data
				istringstream iss(line);
				vector<string> tokens;
				copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter(tokens));

				pair<double, double> tmp_water;

				tmp_water.first = atof(tokens[0].c_str());
				tmp_water.second = atof(tokens[1].c_str());

				this->waterCoords.push_back(tmp_water);
			}
			else if (velocity)
			{
				istringstream iss(line);
				vector<string> tokens;
				copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter(tokens));

				init_velocity.second = atof(tokens[0].c_str());
				init_velocity.first = atof(tokens[1].c_str());
			}
			else if (is_rho) rho_in = atof(line.c_str());

			if (line.find("wet domain:", 0) != string::npos)
				water = true;
			else if (line.find("velocity:", 0) != string::npos)
				velocity = true;
			else if (line.find("density[kg/m3]:", 0) != string::npos)
				is_rho = true;
		}
	}
	else
	{
		assert(bound.is_open());
		printf("file 'BC' is not available");
		assert(init.is_open());
		printf("file 'IC' is not available");
		assert(system.is_open());
		printf("file 'system' is not available");
	}
	assert(waterCoords.size() % 2 == 0 && "please specify a min and max for every square filled with water");

	this->h = this->dx * this->h_fac;
	initialise_grid();

	place_points(init_velocity, rho_in);


	remove_overlap_particles();
	
	//cout << particle_list.size() << endl;	
}


// void SPH_main::set_values(void)
// {
//     min_x[0] = 0.0;
//     min_x[1] = 0.0;

//     max_x[0] = 1.0;
//     max_x[1] = 1.0;

//     h_fac = 1.3;
//     h = dx*h_fac;
// }


void SPH_main::place_points(pair<double, double> init_velocity, double rho_in)
{
	// push the last coord in front and the first coord to the end pof the deque
	this->boundaryCoords.push_front(this->boundaryCoords.back());
	this->boundaryCoords.push_back(this->boundaryCoords[1]);
	// define som deques needed
	deque<double> angles;
	deque<double> angular_lenghths;
	deque<double> dqs;
	deque<double> x_grads;
	deque<double> y_grads;

	for (unsigned int i = 1; i < boundaryCoords.size() - 1; i++)
	{
		// the the vectors that connect the two coords
		double x1 = (boundaryCoords[i].first - boundaryCoords[i-1].first);
 		double y1 = (boundaryCoords[i].second - boundaryCoords[i-1].second);
 		double x2 = (boundaryCoords[i].first - boundaryCoords[i+1].first);
 		double y2 = (boundaryCoords[i].second - boundaryCoords[i+1].second);

		// get the angle between them two vectors
		double angle = acos((x1 * x2 + y1 * y2) / (sqrt(pow(x1, 2) + pow(y1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2))));

		// initiate vars to store the gradiant of the x and y normal at corners 
		double x_grad;
		double y_grad;

		// get the gradient
		if (angle == M_PI/2)
		{
			// special case for 90 deg corners
			x_grad = (x1 + x2) / abs((x1 + x2));
			y_grad = (y1 + y2) / abs((y1 + y2));
		}
		else
		{
			// other cases 
			if (x1/abs(x1) == x2/abs(x2))
				x_grad = x1/abs(x1);
			else
				x_grad = (max(x1/(y1 + 0.0001), x2/(y2 + 0.0001))/abs(max(x1/(y1 + 0.0001), x2/(y2 + 0.0001))));

			if (y1/abs(y1) == y2/abs(y2))
				y_grad = y1/abs(y1);
			else
				y_grad = -(max(x1/(y1 + 0.0001), x2/(y2 + 0.0001))/abs(max(x1/(y1 + 0.0001), x2/(y2 + 0.0001))));
		}
		
		// get the length neede to deviate out the boundaries at the corners 
		double angular_length = 2 * h / cos(M_PI/2 - angle/2);
		// the equivalent of dx in direction if the corners
		double dq = angular_length * dx / (2 * h);

		// push everything to the deques specified above
		angles.push_back(angle);
		angular_lenghths.push_back(angular_length);
		dqs.push_back(dq);
		x_grads.push_back(x_grad);
		y_grads.push_back(y_grad);
	}
	// again push the first elemtn to the back and 
	// the last element to the front if the dque
	angles.push_front(angles.back());
	angles.push_back(angles[1]);

	angular_lenghths.push_front(angular_lenghths.back());
	angular_lenghths.push_back(angular_lenghths[1]);

	dqs.push_front(dqs.back());
	dqs.push_back(dqs[1]);

	x_grads.push_front(x_grads.back());
	x_grads.push_back(x_grads[1]);

	y_grads.push_front(y_grads.back());
	y_grads.push_back(y_grads[1]);
	double cnt = 0;
	// start looping through coords
	for (unsigned int i = 1; i < boundaryCoords.size() - 1; i++)
	{	
		// these define the thinknes of boundary
		double dq1 = 0;
		double dq2 = 0;
		// ´compute the normal vectoir towards of each boundary in direction ti the centre of domain
		double normal2 = boundaryCoords[i+1].first - boundaryCoords[i].first;
		double normal1 = boundaryCoords[i].second - boundaryCoords[i+1].second;
		// normalize the normal vec
		double size = sqrt(normal1 * normal1 + normal2 * normal2);

		normal1 /= size;
		normal2 /= size;

		// loop ober the boundary thikness
		while (dq1 <= angular_lenghths[i])
		{
			//get posotions of start of the boundary at dq
			double x_pos1 = boundaryCoords[i].first + (x_grads[i] * dq1 * cos(angles[i]/2));
			
			double x_pos2 = boundaryCoords[i+1].first + (x_grads[i+1] * dq2 * cos(angles[i+1]/2));

			double y_pos1 = boundaryCoords[i].second + (y_grads[i] * dq1 * sin(angles[i]/2));
			
			double y_pos2 = boundaryCoords[i+1].second + (y_grads[i+1] * dq2 * sin(angles[i+1]/2));

			// get the diference between start and end of upcoming loop
			double diff = sqrt(pow(x_pos1 - x_pos2, 2) + pow(y_pos1 -y_pos2, 2));

			int n = (diff / dx) + 0.9999;
			// loop over lengths of boundary 
			for (int k = 0; k <= n; k++)
			{
				// iniotiate a particle object at every new dx
				SPH_particle tmp;
				// give it some properties
				tmp.x[0] = x_pos1 + (x_pos2 - x_pos1) / n * k;
				tmp.x[1] = y_pos1 + (y_pos2 - y_pos1) / n * k;
				tmp.v[0] = 0;
				tmp.v[1] = 0;
				tmp.P = 0;
				tmp.mass = rho_in * dx * dx * 4;
				tmp.rho = rho_in;
				tmp.norm[0] = normal1;
				tmp.norm[1] = normal2;
				tmp.boundary_particle = true;
				tmp.calc_index();
				this->particle_list.push_back(tmp);
				cnt++;
			}

			dq1 += dqs[i];
			dq2 += dqs[i+1];
		}
	}
	
	// define whether to move right or left from strat point
	double dx_water = (waterCoords[1].first - waterCoords[0].first) /\
	abs(waterCoords[0].first - waterCoords[1].first) * dx;
	// allocate to grid previously computed thingies to grid 
	allocate_to_grid();
	// loop through water coord pairs
	for (int i = 0; i < waterCoords.size() / 2; i++)
	{
		
		
		
		// loop over horizontal length of wet square
		double j = waterCoords[2 * i].first;
		while (j < waterCoords[2 * i + 1].first)
		{
			// loop over vertical length of wet square
			double k = waterCoords[2 * i].second;
			while ( k >= waterCoords[2 * i + 1].second)
			{
				bool to_close_v = false;
				SPH_particle tmp;
				tmp.x[0] = j;
				tmp.x[1] = k;
				tmp.calc_index();
				neighbour_iterate(&tmp);
				for (int n = 0; n < tmp.neighbours.size(); n++)
				{
					double dist = sqrt(pow(tmp.neighbours[n]->x[0] - tmp.x[0], 2) + \
					pow(tmp.neighbours[n]->x[1] - tmp.x[1], 2));
					if (dist < this->dx/2)
						to_close_v = true;

				}

				if (to_close_v)
					break;

				// add deviattion to 
				double deviation = ((double)(rand() % 2000) / 1000 - 1) * dx / 10;
				tmp.x[0] += deviation;

				deviation = ((double)(rand() % 2000) / 1000 - 1) * dx / 10;
				tmp.x[1] += deviation;

				tmp.v[0] = init_velocity.first;
				tmp.v[1] = init_velocity.second;
				tmp.rho = rho_in;
				tmp.P = tmp.cal_pressure(rho_in, this->gamma, this->c0);
				tmp.boundary_particle = false;
				this->particle_list.push_back(tmp);
				cnt++;
				k -= dx_water;
			}
			j += dx_water;
		}
	}
	cout << "cnt: " << cnt << endl;


	// for (int i = 0; i < waterCoords.size() / 2; i++)
	// 	for (double j = waterCoords[2 * i].first; j < waterCoords[2 * i + 1].first; j += dx)
	// 		for (double k = waterCoords[2 * i].second; k < waterCoords[2 * i + 1].second; k += dx)
	// 		{
	// 			SPH_particle tmp;

	// 			double deviation = ((double)(rand() % 2000) / 1000 - 1) * dx / 10;
	// 			tmp.x[0] = j + deviation;

	// 			deviation = ((double)(rand() % 2000) / 1000 - 1) * dx / 10;
	// 			tmp.x[1] = k + deviation;

	// 			tmp.v[0] = init_velocity.first;
	// 			tmp.v[1] = init_velocity.second;
	// 			tmp.rho = rho_in;
	// 			tmp.P = tmp.cal_pressure(rho_in, this->gamma, this->c0);
	// 			tmp.boundary_particle = false;
	// 			tmp.calc_index();
	// 			this->particle_list.push_back(tmp);
 	// 		}
	
}

double SPH_particle::w_function(double dist) {
	//const double M_PI = 3.1415;
	double coef = 10. / (7 * M_PI * pow(SPH_main::h, 2));

	double q = dist / SPH_main::h;

	if (q <= 1) {
		return coef * (1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3));
	}
	else if (q <= 2) {
		return coef * 0.25 * pow(2 - q, 3);
	}
	else {
		return 0;
	}
}

void SPH_particle::update_rho() {
	/* Function that updates a particle's value of rho.
	 This function sould be called every 10 or 20 steps in the iteration/ It also needs to be run inside a loop over every particle when called.  */

	 // the calculation of updating the density value uses the cubic spline function (W(r,h)) and the density of the neighbour

	double top_calc = 0.0; // top calculation of equation
	double bottom_calc = 0.0; // bottom calculation of equation
	double dn[2];
	double dist;

	// it needs to know all of its neighbours, and loop through every neighbour
	for (int i = list_num[0] - 1; i <= list_num[0] + 1; i++) 
		if (i >= 0 && i < main_data->max_list[0]) 
			for (int j = list_num[1] - 1; j <= list_num[1] + 1; j++) 
				if (j >= 0 && j < main_data->max_list[1])
				{
					//#pragma omp parallel for
					for (unsigned int cnt = 0; cnt < main_data->search_grid[i][j].size(); cnt++)
					{
						SPH_particle *other_part = main_data->search_grid[i][j][cnt];
						if (this != other_part)                    //stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
								dn[n] = x[n] - other_part->x[n];
							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

							if (dist < 2. * main_data->h)                    //only particle within 2h
							{
								top_calc += w_function(dist);
								bottom_calc += (w_function(dist) / other_part->rho);
							}
						}
					}
				}
	this->rho_new = top_calc / bottom_calc;
}


void SPH_main::initialise_grid(void)
{
	min_x[0] = boundaryCoords[0].first;
    min_x[1] = boundaryCoords[0].second;

    max_x[0] = boundaryCoords[0].first;
    max_x[1] = boundaryCoords[0].second;
	
	for (int i = 1; i < boundaryCoords.size(); i++)
	{
		if (boundaryCoords[i].first < min_x[0])
			min_x[0] = boundaryCoords[i].first;
		if (boundaryCoords[i].second < min_x[1])
			min_x[1] = boundaryCoords[i].second;
		if (boundaryCoords[i].first > max_x[0])
			max_x[0] = boundaryCoords[i].first;
		if (boundaryCoords[i].second > max_x[1])
			max_x[1] = boundaryCoords[i].second;
	}

	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 2.0 * h;
		max_x[i] += 2.0 * h;                                              //add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1);
	}
	search_grid.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
		search_grid[i].resize(max_list[1]);
}

void SPH_particle::calc_index()
{
	for (int i = 0; i < 2; i++)
	{
		list_num[i] = int((this->x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
	}
}

void SPH_main::allocate_to_grid(void)                //needs to be called each time that all the particles have their positions updated
{
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
	
}


void SPH_main::neighbour_iterate(SPH_particle* part)                    //iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
	SPH_particle* other_part;
	double dist;            //distance between particles
	double dn[2];            //vector from 1st to 2nd particle

	double C = 20;

	part->a[0] = 0.0;
	part->a[1] = -9.81;
	

	part->D = 0;

	// double refl_x = 0;
	// double refl_y = 0;
	bool dir_neigh = false;

	part->neighbours.clear();
//	#pragma omp parallel for schedule(static)
	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					//#pragma omp parallel for
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];
						if (part != other_part)                    //stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
								dn[n] = part->x[n] - other_part->x[n];
							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

							if (dist < dx && other_part->boundary_particle && !part->boundary_particle)
								dir_neigh = true;

							// if (other_part->boundary_particle && !part->boundary_particle)
							// {
							// 	//#pragma omp atomic
							// 	refl_x += dn[0];
							// 	//#pragma omp atomicƒ
							// 	refl_y += dn[1];
							// }

							if (dist < 2. * h)                    //only particle within 2h
							{
								//TODO: all the interactions between the particles
								part->neighbours.push_back(other_part);
								part->navier(other_part, dx, gamma, c0, mu);

								//cout << "dn: " << dn[0] << " " << dn[1] << endl;        //Should be removed from the code - simply here for you to see that it is working
							}
						}
					}
				}
	// part->a[0] += pow(refl_x * part->main_data->dx , gamma) / C;
	// part->a[1] += pow(refl_y * part->main_data->dx , gamma) / C;
	//cout << refl_y << " "<< 3 * h/dx<<endl;

	// if (dir_neigh)
	// {
	// 	part->x[0] += other_part->norm[0] * dx;
	// 	part->x[1] += other_part->norm[0] * dx;
	// }

	if (dir_neigh)
	{
		part->x[0] += -part->v[0] * other_part->norm[0] * dt;
		part->v[0] = -part->v[0];
	}

	if (dir_neigh)
	{
		part->x[1] += -part->v[1] * other_part->norm[0] * dt;
		part->v[1] = -part->v[1];

	}
	
}

double SPH_main::dW_dr(double r)
{
	double q = r / h;

	// if (q >= 0 && q <= 1)
	//     return (10 / (7 * M_PI * pow(h, 2))) * (1 - 3/2 * pow(q, 2) + 3/4 * pow(q, 3));
	// else if (q >= 1 && q <= 2)
	//     return (10 / (7 * M_PI * pow(h, 2))) * (1/4 * pow(2 - q, 3));
	// else //if (q > 2)
	//     return 0; 
	//const double M_PI = 3.1415;
	if (q <= 1)
		return (10. / (7. * M_PI * pow(h, 3))) * (-3 * q + (9. / 4.) * pow(q, 2.0));
	else if (q <= 2)
		return -(10. / (7. * M_PI * pow(h, 3))) * (0.75 * pow(2 - q, 2.0));
	else
		return 0;
}

void SPH_main::remove_overlap_particles(/*vector<SPH_particle> particle_list*/)
{
	/* a function which checks for concurrent or near concurrent points and removes one of them */
// first loop through every particle
//double initial_num_particles = particle_list.size();

	for (int i = 0; i < particle_list.size(); i++) // less than or less than and equal to?
	{
		// have to find all the neighbours
		neighbour_iterate(&particle_list[i]);

		// loop through particle's neighbours
		double dn[2];            //vector from 1st to 2nd particle
		for (int n = 0; n < particle_list[i].neighbours.size(); n++)
		{
			for (int l = 0; l < 2; l++)
			{
				dn[l] = particle_list[i].neighbours[n]->x[l] - particle_list[i].x[l];
			}
			double dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

			// find the distance between the current particle and current neighbour

		// if the distance between the current particle and the neighbour is zero or less than a certain value (i.e. close to overlap) then delete one of the particles.
			// WHAT should this value be??? 1% of initial particle spacing dx??
			if (dist <= dx / 2)
			{
				if (particle_list[i].boundary_particle == false) // making sure i dont delete a particle that is a boundary particle
				{
					// delete the particle
					particle_list.erase(particle_list.begin() + i);
					// does this affect the outer loop where i < particle_list.size() because now the size is smaller?
					// fixed this by finding the initial number of particles
					//
					//i--; // ? because if you delete 2nd element in list of particles for example, the particle after it will now become the second element in the list? so it will accidentally get skipped in the loop - so minusing one from list so it isnt skipped over
				}
			}
		}
	}
}
