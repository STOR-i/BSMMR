#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include "modelfunctions.hpp"
#include "calculate_accept_prob.hpp"
#include "function_evaluation.hpp"

using namespace std;

#ifndef UPDATES_HPP
#define UPDATES_HPP

/*Sample the proposal, calculate the acceptance probability and update process accordingly*/

void update_process(const int &k, const int &ip, 
		    vector<vector<vector<vector<double> > > > &processes,
		    const vector<vector<int> > &ind, const vector<vector<double> > &covariates_k,
		    vector<double> &assigned_level_current, const double &baseline_k, 
		    const vector<double> &response_k, const vector<int> &number_trials_k,
		    const double &param_k, const string &set, const vector<int> &neighbours,
		    const vector<double> &weights,
		    const vector<vector<double> > & sample_boundaries,
		    const vector<vector<vector<double> > > &A_k,
		    const vector<vector<int> > &low_sets, const vector<vector<int> > &upp_sets,
		    const vector<vector<vector<int> > > &pos_cov,
		    const double &volume_sample_space,
		    const double &DELTA_MIN, const double &DELTA_MAX,
		    const double &omega, const int &dim, const double &p, const double &q,
		    const double &eta, const double &p_birth, const double &p_death,
		     double &sum_levels, int &number_jumps, int &number_jumps_subprocess,
		    const int & maximum_jumps, mt19937 &generator,
		    void (*model_log_like)(const double &, const double &, const double &,
					   const double &, const int &, const double &, double &));

/*Update the random effects for all K regions*/

void update_baseline(vector<double> &baseline, const int&K,const vector<vector<double> > &response,
		     const vector<vector<int> > &number_trials,
		     const vector<vector<double> > &assigned_level,
		     const vector<vector<int> > &neighbours,const vector<vector<double> > &weights,
		     const vector<double> &sum_weights, vector<double> &mean,
		     const vector<double> &param, const double &tau, mt19937 &generator,
		     void (*model_log_like)(const int &, const vector<double> &, const double &,
					    const double &, const vector<double> &,
					    const vector<int> &, const double &, double &));

/*Update the spatial smoothing parameter for the random effects*/

void update_tau(double &tau, const vector<double> &sum_weights,
		const vector<double> &random_effects, const vector<vector<double> > weights,
		const vector<vector<int> > neighbours, const int &K,
		const double &alpha_tau, const double &beta_tau, mt19937 &generator);

/*Update additional parameter in the distribution*/

void update_param(const string &model, vector<double> &param, const int &K,
		  const vector<vector<double> > &level, const vector<double> &random_effect,
		  const vector<vector<double> > &number_claims, const double &alpha_param,
		  const double &beta_param, mt19937 &generator);

/*Sample proposed support point and associated level for a subprocess*/

vector<vector<double> > sample_proposal_birth(const int &ip,
					      const vector<vector<vector<double> > > &process,
					      const vector<vector<int> > &ind,
					      const vector<vector<double> >&sample_boundaries,
					      const double &DELTA_MIN, const double &DELTA_MAX,
					      const vector<int> &low, const vector<int> &upp,
					      const vector<vector<vector<int > > > &pos_cov,
					      mt19937 &generator);

/*Sample support point and associated level to be proposed for being removed from the subprocess*/


vector<vector<double> > sample_proposal_death(const int &ip,
					      vector<vector<vector<double> > > &process,
					      const vector<vector<int> > &ind,
					      const double &DELTA_MIN, const double &DELTA_MAX,
					      const vector<int> &low, const vector<int> &upp,
					      const vector<vector<vector<int > > > &pos_cov, 
					      mt19937 &generator);

/*Sample proposed support point and associated level for a subprocess for a shift being proposed*/


vector<vector<double> > sample_proposal_shift(const int &ip,
					      vector<vector<vector<double> > > &process,
					      const vector<vector<int> > &ind,
					      const vector<vector<double> > &sample_boundaries,
					      const double &DELTA_MIN, const double &DELTA_MAX,
					      const vector<int> &low, const vector<int> &upp,
					      const vector<vector<vector<int > > > &pos_cov,
					      double &sum_levels, mt19937 &generator);

#endif
