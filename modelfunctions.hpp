#include <iostream>

using namespace std;

#ifndef MODEL_FUNCTIONS_HPP
#define MODEL_FUNCTIONS_HPP

/*Likelihood ratio on log scale for Binomial model with success probability on logit scale
  modelled by monotonic function*/

void Binomial_model(const double &level_proposed, const double &level_current,
		    const double &baseline, const double &response, const int &number_trials,
		    const double &param, double &result);

/*Likelihood ratio on log scale for Poisson model with rate on log scale modelled by monotonic
  function*/

void Poisson_model(const double &level_proposed, const double &level_current,
		   const double &baseline, const double &response, const int &number_trials,
		   const double &param, double &result);

/*Likelihood ratio on log scale for Gaussian model with mean modelled by monotonic function*/

void Gaussian_model(const double &level_proposed, const double & level_current,
		    const double &baseline, const double &response, const int &number_trials,
		    const double &variance, double &result);

/*Likelihood ratio on log scale for Binomial model truncated at 0 with success probability on
  logit scale modelled by monotonic function*/

void trBinomial_model(const double &level_proposed, const double &level_current,
		      const double &baseline, const double &response, const int &number_trials,
		      const double &param, double &result);

/*Likelihood ratio on log scale for Generalized Pareto model with scale on log scale modelled
  by monotonic function*/

void GPD_model(const double &level_proposed, const double &level_current,
	       const double &baseline, const double &response, const int &number_trials,
	       const double &shape, double &result);

/*Likelihood ratio on log scale for Bernoulli model with succes probability on logit scale modelled
  by monotonic function*/

void Bernoulli_model(const double &level_proposed, const double &level_current,
		     const double &baseline, const double &response,
		     const int &number_trials, const double &param, double &result);

/*Likelihood ratio on log scale for Binomial model with success probability on logit scale
 modelled by monotonic function for proposed new baseline*/

void Binomial_model_re(const int &n, const vector<double> &level, const double &baseline_current,
		       const double &baseline_proposed, const vector<double> &response,
		       const vector<int> &number_trials, const double &param, double &result);

/*Likelihood ratio on log scale for Poisson model with rate on log scale modelled by monotonic
 function for proposed new baseline*/

void Poisson_model_re(const int &n, const vector<double> &level, const double &baseline_current,
		      const double &baseline_proposed, const vector<double> &response,
		      const vector<int> &number_trials, const double &param, double &result);

/*Likelihood ratio on log scale for Gaussian model with mean modelled by monotonic function
 for proposed new baseline*/

void Gaussian_model_re(const int &n, const vector<double> &level, const double &baseline_current,
		       const double &baseline_proposed, const vector<double> &response,
		       const vector<int> &number_trials, const double &variance, double &result);

/*Likelihood ratio on log scale for Binomial model truncated at 0 with success probability on
  logit scale modelled by monotonic function for proposed new baseline*/

void trBinomial_model_re(const int &n, const vector<double> &level, const double &baseline_current,
			 const double &baseline_proposed, const vector<double> &response,
			 const vector<int> &number_trials, const double &param, double &result);

/*Likelihood ratio on log scale for Generalized Pareto model with scale on log scale modelled
  by monotonic function for proposed new baseline*/

void GPD_model_re(const int &n, const vector<double> &level, const double &baseline_current,
		  const double &baseline_proposed, const vector<double> &response,
		  const vector<int> &number_trials, const double &shape, double &result);

/*Likelihood ratio on log scale for Bernoulli model with succes probability on logit scale modelled
  by monotonic function for proposed new baseline*/


void Bernoulli_model_re(const int &n, const vector<double> &level, const double &baseline_current,
			const double &baseline_proposed, const vector<double> &response,
			const vector<int> &number_trials, const double &param, double &result);

#endif
