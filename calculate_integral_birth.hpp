#include <iostream>
#include <vector>

#ifndef CHANGE_INTEGRAL_BIRTH
#define CHANGE_INTEGRAL_BIRTH

using namespace std;

/*Calculate prior ratio for proposing a 'Birth' for Region k*/

double integral_birth(const string &set, const int &k, const vector<double> &proposal,
		      const int &ip, const vector<vector<vector<vector<double> > > > &processes,
		      const vector<vector<int> > &ind, const vector<int> &neighbours,
		      const vector<double> &weights, const int &dim, const double &DELTA_MIN,
		      const vector<vector<double> > &sample_boundaries,
		      const vector<vector<vector<double> > > &A_k, const vector<int> &low_sets,
		      const vector<vector<vector<int> > > &pos_cov, const double &level_curr,
		      const double &p, const double &q);

/*Derive corner points of the area affected by the proposal and the current function over it*/

void simplify_process(vector<double> &proposal_full, const vector<double> &proposal,
		      const vector<vector<vector<double> > > &process,
		      const vector<vector<int> > &ind, const int &dim, const double &level_curr, 
		      const vector<vector<int> > &pos_cov, vector<vector<double> > &process_k,
		      vector<vector<double> > &border, vector<vector<double> > &temp_border,
		      vector<double> &bound);
/*Derive set of points describing area affected by the proposal from set of potential candidates*/

void minimum_border(vector<vector<double> > &border, vector<vector<double> > &temp_border,
		    vector<double> &bound, const int &dim);

/*Derive set of points describing current function over the area affected by the proposal*/

void minimum_process_check_border(vector<vector<double> > &process, 
				  const vector<vector<double> > &border, const int &dim);

/*Delete points which are not describing the current function over the area affected the proposal*/

void minimum_process(vector<vector<double> > &process, const int &dim, const int &pos);

/*Check whether a candidate descrbining the current function over the area affected by the
  proposal lies outside the area by checking monotonicity for all but the first pos border points*/

bool check_border( const vector<double> &point, const vector<vector<double> > &border,
		   const int &dim, const int &pos);

/*Sort points with respect to their mark*/

void sort_by_level(vector<vector<double> > &process);

/*Derive function over an area given the border points and the proposal for neighbour*/

void derive_process_neigh(vector<double> &proposal_full, const vector<double> &proposal,
			  const vector<vector<vector<double> > > &process,
			  const vector<vector<int> > &ind, const int &dim,
			  const double &temp_level, const vector<vector<int> > &pos_cov,
			  const vector<vector<double> > &border, const vector<double> &bound,
			  const vector<vector<double> > &temp_border, 
			  vector<vector<double> > &process_k2);

/*Calculate the integrated difference of the process the proposal and its neighbours*/

void integral_difference(vector<vector<double> > &process_k, 
			 vector<vector<vector<double> > > &processes_neigh,
			 vector<vector<double> > &border, vector<vector<double> > &temp_border_old,
			 const int &dim, const double &proposal, vector<int> ind, 
			 const double &DELTA_MIN, const double &p, const double &q,
			 const vector<double> &weights, double &integral);

/*Derive set of border points describing the area with level equal to level at a point*/

void derive_temp_border(const vector<double> &point,vector<vector<double> > &border,
			const int &dim, vector<int> &ind, vector<double> &proposal,
			vector<vector<double> > &temp_border);

/*Derive maximum set of points such that no point puts monotonic constraint on another point*/

void minimum_region(vector<vector<double> > &border, const int &dim);

/* Calculation of the volume of a multi-dimensional polytope with given border points*/

double calculate_volume(vector<vector<double> > &process_border, int dim);

#endif
