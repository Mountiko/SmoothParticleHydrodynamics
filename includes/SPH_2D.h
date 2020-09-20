#pragma once
#include <vector>
#include <deque>
#include <cmath>
#include <iostream>
#include <utility>

using namespace std;

class SPH_main;

class SPH_particle
{
public:

    bool boundary_particle;
    double x[2], v[2];                //position and velocity
    double a[2];                     //acceleration
    double P = 0;                    //density and pressure
    double mass;
    double rho, rho_new;	//density
    double D;
    vector<SPH_particle*> neighbours;

	//double norm[2];

    static SPH_main* main_data;        //link to SPH_main class so that it can be used in calc_index

    int list_num[2];                //index in neighbour finding array

    void calc_index();
    double cal_mass(double dx, double rho);

	double norm[2];

    double cal_pressure(double rho, double gamma, double c0);

    //function to calculate dv/dt
    void navier(SPH_particle* neigh, double dx, double gamma, double c0, double mu);
    //function to calculate drho/dt
    double w_function(double dist);
    void update_rho(void);
};


class SPH_main
{
public:

    static double h;                                //smoothing length

    SPH_main(); // constructor
	SPH_main(SPH_main* domain);
    void initialise_grid(void);


    void set_values(void);

    void place_points(pair<double, double> init_velocity, double rho_in);

    void allocate_to_grid(void);            //allocates all the points to the search grid (assumes that index has been appropriately updated)

    void neighbour_iterate(SPH_particle* part);

    void remove_overlap_particles(/*vector<SPH_particle> particle_list*/);

	void fill_water(double x1, double x2);

    double c0;
    double gamma;

    double mu; // viscosity

    double sim_time; // time of simulation

    static double dW_dr(double r);

	double dt;

    double h_fac;
    double dx;                                //particle initial spacing

    double min_x[2], max_x[2];                //dimensions of simulation region

	string time_step_meth;

	bool adaptive_timestepping;

    int max_list[2];

    vector<pair<double, double>> waterCoords;
    deque<pair<double, double>> boundaryCoords;

    vector<SPH_particle> particle_list;				//list of all the particles

    vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell
};

