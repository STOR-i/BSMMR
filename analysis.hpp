#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <random>
#include "function_evaluation.hpp"

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

/*Read MCMC output and call further functions to analyse the MCMC output and predict test data*/

void analysis_and_prediction(const string &model, const int&dim, const int &K,
			     const int &grid_size, const double &DELTA_MIN,
			     const int &convergence_samples_size, 
			     const vector<vector<vector<double> > > &sample_boundaries, 
			     mt19937 &generator);

/*Check convergence by sampling a given number of points and derive assigned assigned levels
  for these points for each iteration step*/

void check_convergence(const int &sample_size, const vector<vector<vector<double> > > &process_k,
		       const vector<double> &baseline_k, const int &k, const int &dim,
		       const double &DELTA_MIN, const vector<vector<double> > &sample_boundaries,
		       mt19937 &generator);

/*Derive mean and quantiles of the estimated posterior function on a regular grid for Region k*/ 

void function_on_grid(const vector<vector<vector<double> > > &process,
		      const vector<double> &baseline, const int &dim, const int &k,
		      const int &grid_size, const double &DELTA_MIN,
		      const vector<vector<double> > &sample_boundaries);

/*Derive next grid point to be considered*/

void update_pos(vector<int> &pos, const int &j, const int &grid_size, bool &end_of_grid);

/*Derive mean and quantiles for a given grid point*/

vector<double> derive_summary_statistics(const vector<double> &point, 
					 const vector<vector<vector<double > > > &process_k, 
					 const vector<double> &baseline_k, const int &dim,
					 const double &DELTA_MIN);

/*Detection of potential threshold effects for Region k from the MCMC chains*/

void detect_thresholds(const vector<vector<vector<double> > > &process_k, 
		       const vector<double> &baseline_k, const vector<double> &min_element,
		       const int &dim, const int &k, const double &DELTA_MIN);

/* Import test data set from external "test_data.txt" file*/

void import_test_data(vector<vector<vector<double> > > &test_covariates,
		      vector<vector<int> > &test_number_trials,
		      vector<vector<double> > &true_response,
		      const int &dim, const string &model);

/* Derive mean predictive error for a set of test points*/

void prediction(const vector<vector<vector<double> > > &process, const vector<double> &baseline,
		const vector<double> &param, const vector<vector<double> > &test_covariates,
		const vector<int> &test_number_trials, const vector<double> &true_response,
		const int &dim, const string &model, const double &DELTA_MIN);

/*Derive predictive mean for multiple points based on their assigned levels from the MCMC
  algorithm*/

void derive_predictive_mean(const string &model, const vector<vector<double> > &levels,
			    const vector<double> &param, const vector<int> &number_trials,
			    vector<double> &predictive_mean);

#endif
