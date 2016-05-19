#include <vector>
#include <iostream>

using namespace std;

#ifndef CALCULATE_ACCEPTANCE_PROCESS_HPP
#define CALCULATE_ACCEPTANCE_PROCESS_HPP

/*Calculation of the acceptance probability for a proposed birth*/

double accept_prob_birth(const int &k, vector<vector<vector<vector<double> > > > &processes, 
			 const vector<vector<int> > &ind, const vector<vector<double> > &proposal,
			 const int& ip, const vector<vector<double> > &covariates_k, 
			 const vector<double> &assigned_level_curr,
			 vector<double> &assigned_level_prop, const double &baseline_k,
			 const vector<double> &response_k, const vector<int> &number_trials_k,
			 const double &param_k, const string &set, const vector<int> &neighbours,
			 const vector<double> &weights, 
			 const vector<vector<double> > &sample_boundaries,
			 const vector<vector<vector<double> > > &A_k, const vector<int> &low_sets,
			 const vector<vector<vector<int> > > &pos_cov, const double &DELTA_MIN,
			 const double &omega, const int &dim, const double &p, const double &q,
			 const double &eta,  const double &p_birth, const double &p_death,
			 void (*model_log_like)(const double &, const double &,
						const double &, const double &, const int&,
						const double &, double &));
	      
/*Calculation of the acceptance probability for a proposed death*/

double accept_prob_death(const int &k, vector<vector<vector<vector<double> > > > &processes, 
			 const vector<vector<int> > &ind, const vector<vector<double> > &proposal,
			 const int& ip, const vector<vector<double> > &covariates_k, 
			 const vector<double> &assigned_level_curr,
			 vector<double> &assigned_level_prop, const double &baseline_k,
			 const vector<double> &response_k, const vector<int> &number_trials_k,
			 const double &param_k, const string &set, const vector<int> &neighbours,
			 const vector<double> &weights, 
			 const vector<vector<double> > &sample_boundaries,
			 const vector<vector<vector<double> > > &A_k, const vector<int> &low_sets,
			 const vector<vector<vector<int> > > &pos_cov, const double &DELTA_MIN,
			 const double &omega, const int &dim, const double &p, const double &q,
			 const double &eta,  const double &p_birth, const double &p_death,
			 void (*model_log_like)(const double &, const double &,
						const double &, const double &,
						const int &, const double &, double &));
   
/*Calculation of the acceptance probability for a proposed "Shift"*/

double accept_prob_shift(const int &k, vector<vector<vector<vector<double> > > > &processes, 
			 const vector<vector<int> > &ind, const vector<vector<double> > &proposal,
			 const int &ip, const vector<vector<double> > &covariates_k, 
			 const vector<double> &assigned_level_curr,
			 vector<double> &assigned_level_prop, const double &baseline_k,
			 const vector<double> &response_k, const vector<int> &number_trials_k,
			 const double &param_k, const string &set, const vector<int> &neighbours,
			 const vector<double> &weights,
			 const vector<vector<double> > &sample_boundaries,
			 const vector<vector<vector<double> > > &A_k, const vector<int> &low_sets,
			 const vector<vector<vector<int> > > &pos_cov, const double &DELTA_MIN,
			 const double &omega, const int &dim, const double &p,
			 const double &q,const double &sum_levels, const double &temp_sum_levels,
			 void (*model_log_like)(const double &, const double &,
						const double &,	const double &, const int&,
						const double &, double &));

/*Derive likelihood ratio on log-scale for a proposed "Shift" and update levels for data points*/

void calculate_likelihood_shift(double &result, const vector<vector<double> > &proposal,
				const int &ip, const vector<vector<vector<double> > > &process,
				const vector<vector<int> > &ind, const int &dim,
				const double &baseline_k,
				const vector<vector<double> > &covariates_k,
				const vector<double> &assigned_level_curr, 
				vector<double> &assigned_level_prop,
				const vector<double> &response_k,
				const vector<int> &number_trials_k, const double &param_k,
				void(*model_log_like)(const double &, const double &,
						      const double &, const double &, const int &,
						      const double &, double &));

/*Calculation of the acceptance probability for a proposed baseline*/

double accept_prob_baseline(const vector<double> &assigned_level_curr,
			    const vector<double> &response_k, const double &baseline_prop,
			    const double &baseline_curr, const vector<int> &number_trials_k,
			    const double &param_k, const double &mean, const double &sum_weights,
			    const double &rho,
			    void (*model_log_like)(const int &, const vector<double> &,
						   const double &, const double &,
						   const vector<double> &, const vector<int> &,
						   const double &, double &));

#endif
