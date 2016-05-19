#include <vector>
#include <iostream>
#include "function_evaluation.hpp"

#ifndef UTILITIES_MONOTONIC_REGRESSION_HPP
#define UTILITIES_MONOTONIC_REGRESSION_HPP

/*Read the parameters for RJMCMC algorithm, model and analysis*/

void read_configurationfile(int &number_iterations, int &burn_in, int &size_thinning,
			    double &p_birth, double &p_death, int &max_jumps, string &model,
			    double &alpha_param, double &beta_param, int &K, int &dim,
			    double &DELTA_MIN, double &DELTA_MAX, string &set, string &cross_val,
			    double &omega, double &omega_0, int &number_repetitions,
			    double &p, double &q, double &eta, double &tau_0, double &alpha_tau,
			    double &beta_tau, int &grid_size, int &convergence_samples_size);

/*Import training data to fit the model*/

void import_data(vector<vector<vector<double> > > &covariates, vector<vector<double> > &response,
		 vector<vector<int> > &number_trials, vector<vector<int> > &neighbours,
		 vector<vector<double> > &assigned_level_curr, vector<vector<double> > &weights,
		 vector<vector<double> > &weights_baseline, const int &dim,
		 const double &DELTA_MIN, const string &model, const int &K);

/*Initialize the parameters to run the RJMCMC algorithm*/

void initialise_parameters(vector<vector<vector<vector<double> > > > &processes,
			   vector<vector<int> > &ind, vector<double> &baseline,
			   vector<double> &param, vector<vector<int> > &low_sets,
			   vector<vector<int> > &upp_sets, vector<vector<vector<int> > > &pos_cov,
			   vector<vector<double> > &sum_levels,
			   vector<vector<int> > &number_jumps_subprocesses,
			   vector<int> &number_jumps,vector<double> &mean,
			   vector<double> &sum_weights_baseline,  
			   const vector<vector<double> > &weights_baseline, const int &dim,
			   const int &K);

/*Derive lower left and upper right corner of the spaces on which the proposed new points are
  sampled and derive the associated voume of the sample space in dependency on the considered
  domain of integration*/

void derive_volume(vector<vector<vector<vector<double> > > > &A,
		   vector<vector<vector<double> > > &sample_boundaries,
		   vector<vector<double> > &volume_sample_space,
		   const vector<vector<vector<double> > > &covariates,
		   const int &dim, const int &K, const string &set,
		   const vector<vector<int> > &neighbours, const vector<vector<int> > &ind);

/*Store the current processes, random effects and smoothing parameter to different files*/

void write_to_file(const vector<vector<vector<vector<double> > > > &processes, const int &K,
		   const int &dim, const vector<vector<int> > &ind,
		   const vector<double> &baseline, const double &tau, const vector<double> &param,
		   const string &model, const vector<int> &number_jumps,
		   const vector<vector<vector<double> > > &sample_boundaries);

#endif
