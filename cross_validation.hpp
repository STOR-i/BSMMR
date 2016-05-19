#include <vector>
#include <iostream>
#include <random>
#include "updates.hpp"

using namespace std;

#ifndef CROSS_VALIDATION_HPP
#define CROSS_VALIDATION_HPP

/*Derive appropriate value for omega by efficient global optimization and m-fold cross validation*/

double optimal_omega(const double &omega_0, const int &rep, const int &number_folds, const int &K, 
		     const string &model, const vector<vector<int> > &neighbours,
		     const string &set, const vector<vector<double> > &weights,
		     const vector<vector<vector<double> > > &sample_boundaries,
		     const vector<vector<vector<vector<double> > > > A,
		     const vector<vector<double> > &volume_sample_space, const double &DELTA_MIN, 
		     const double &DELTA_MAX, const double &p, const double &q, const int &dim, 
		     const double &eta, const double &tau, const double &alpha_tau,
		     const double &beta_tau, const double &alpha_param, const double &beta_param,
		     const double &p_birth, const double &p_death, const int &max_jumps,
		     const vector<vector<vector<double> > > &covariates,
		     const vector<vector<int> > &number_trials,
		     const vector<vector<double> > &response,
		     const vector<vector<double> >&weights_baseline,
		     void (*model_log_like)(const double &, const double &, const double &,
					    const double &, const int &, const double &, double &),
		     void (*model_log_like2)(const int &, const vector<double> &, const double &,
					     const double &, const vector<double> &,
 					     const vector<int> &, const double &, double &));

/*Perform m-fold cross validation for a given omega and a number of repetitions*/

double cross_validation(const double &omega, const int &rep, const int &number_folds, const int &K,
	    		const string &model, const vector<vector<int> > &neighbours,
			const string &set, const vector<vector<double> > &weights,
			const vector<vector<vector<double> > > &sample_boundaries,
			const vector<vector<vector<vector<double> > > > & A,
			const vector<vector<double> > &volume_sample_space,
			const double &DELTA_MIN, const double &DELTA_MAX, const double &p,
			const double &q, const int &dim, const double &eta, const double &tau,
			const double &alpha_tau, const double &beta_tau, const double &alpha_param,
			const double &beta_param, const double &p_birth, const double &p_death,
			const int &max_jumps, const vector<vector<vector<double> > > &covariates,
			const vector<vector<int> > &number_trials, 
			const vector<vector<double> > &response,
			const vector<vector<double> > &weights_baseline,
			void (*model_log_like)(const double&, const double&, const double&,
					       const double&, const int&, const double&, double&),
			void (*model_log_like2)(const int &, const vector<double>&, const double &,
						const double &, const vector<double> &,
						const vector<int> &, const double &, double &));

/*Split of the original data into training and test data for cross-validation*/

void generate_cv_data(vector<vector<vector<double> > > &test_covariates,
		      vector<vector<int> > &test_number_trials,
		      vector<vector<double> > &test_response,
		      vector<vector<vector<double> > > &training_covariates,
		      vector<vector<int> > & training_number_trials,
		      vector<vector<double> > &training_response,
		      const vector<vector<vector<double> > > &covariates,
		      const vector<vector<double> > &response,
		      const vector<vector<int> > &number_trials, const int &run, const int &fold,
		      const int &K, const int &number_folds, const string &model);

/*Calculate posterior mean for a set of test data points for one region*/

void derive_expectation(const string &model, const vector<vector<double> > &prediction,
			const double &param, const vector<int> &number_trials,
			vector<double> &expectation);

#endif
