#include <vector>
#include <iostream>

#ifndef FUNCEVAL_HPP
#define FUNCEVAL_HPP

using namespace std;

/*Derve functional level of the monotonic function at point arg*/

double function_evaluation(const vector<vector<vector<double> > > &process,
			   const vector<vector<int> > &ind, const int &dim,
			   const vector<double> &arg, const double &level);

/*Derive functional levels for a set of points*/

vector<double> function_evaluation_set(const vector<vector<vector<double> > > &process,
				       const vector<vector<int> > &ind, const int &dim,
				       const vector<vector<double> >&args,
				       const double &DELTA_MIN);

/*Check whether monontonic criterion is fulfilled for all but the first pos components*/

bool check_monotonicity(const vector<double> &point, const vector<double> &arg,
			const vector<int> &ind, const int&pos);

/*Derive the minimum and maximum for sampling a proposed level at sampled location*/

vector<double> bounds_level(const vector<vector<vector<double> > > &process,
			    const vector<double> &proposal, const int &ip,
			    const double &DELTA_MIN, const double &DELTA_MAX,
			    const vector<int> &low_sets, const vector<int> &upp_sets,
			    const vector<vector<vector<int> > > &pos_cov);

/*Check monontonic constraint for two points from potentially different subprocesses*/

bool monotonicity_sets(const vector<double> &point1, const vector<double> &point2,
		       const int &pos, const vector<int>& pos_cov);

/*Derive the current value of the function at the proposal */ 

double level_at_proposal(const vector<vector<vector<double> > > &process,
			 const vector<double> &proposal, const int &ip, const double &DELTA_MIN,
			 const vector<int> &low_sets,const vector<vector<vector<int> > > &pos_cov);

/* Derive sample boundaries for new location for proposed switch. Boundaries are
   derived such that the monotonicity structure is maintained.*/

vector<vector<double> > bounds_location(const int &ip, const vector<double> &point,
					const vector<vector<vector<double> > > &process,
					const vector<vector<int> > &ind, 
					const vector<vector<double> > &sample_boundaries,
					const vector<int> &lows, const vector<int> &upp,
					const vector<vector<vector<int> > > &pos_cov);

/*Derive level at a given point for a process described on the full covariate space only*/

double evaluate_at_point(const vector<double> &point, const vector<vector<double> > &process,
			 const int &dim);

/*Check whether point2 puts constraint on Point 1 with both defined on full covariate space*/

bool monotonicity_full(const vector<double> &point1, const vector<double> &point2,
		       const int &dim);

/*Derive level for a point with marked point process defined on full covariate space*/

double evaluate_test_point(const vector<double> &point, const vector<vector<double> > &process,
			   const double &baseline, const int &dim, const double &DELTA_MIN);

#endif
