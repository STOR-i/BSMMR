#include <iostream>
#include <cmath> 
#include <algorithm>
#include "modelfunctions.hpp"

using namespace std;

/*Likelihood ratio on log scale for Binomial model with success probability on logit scale
 modelled by monotonic function

INPUT:
 level_proposed - Proposed level of the point
 level_current  - Currently assigned level of the point
 baseline       - Currently assigned baseline level
 response       - Number of successes
 number_trials  - Number of trials
 param          - Empty entry
 result         - Current likelihood difference
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void Binomial_model(const double &level_proposed, const double &level_current,
		    const double &baseline, const double &response, const int &number_trials,
		    const double &param, double &result){

  result = result + (level_proposed - level_current) * response - number_trials * 
    ( log(1 + exp(level_proposed + baseline)) - log(1 + exp(level_current + baseline)) );
}

/*Likelihood ratio on log scale for Poisson model with rate on log scale modelled by monotonic
  function

INPUT:
 level_proposed - Proposed level of the point
 level_current  - Currently assigned level of the point
 baseline       - Currently assigned baseline level
 response       - Number of cases
 number_trials  - Empty
 param          - Empty entry
 result         - Current likelihood difference
OUTPUT:
 Updated likelihood ratio of the proposed and current point process
*/

void Poisson_model(const double &level_proposed, const double &level_current,
		   const double &baseline, const double &response, const int &number_trials,
		   const double &param, double &result){
  result = result + response * (level_proposed - level_current) -
    ( exp(level_proposed + baseline) -  exp(level_current + baseline) );
}

/*Likelihood ratio on log scale for Gaussian model with mean modelled by monotonic function

INPUT:
 level_proposed - Proposed level of the point
 level_current  - Currently assigned level of the point
 baseline       - Currently assigned baseline level
 response       - Observation of the response variable
 number_trials  - Empty entry
 variance       - Estimated variance of the Gaussian random variable (constant) 
 result         - Current likelihood difference
OUTPUT:
 Updated likelihood ratio of the proposed and current point process
*/

void Gaussian_model(const double &level_proposed, const double & level_current,
		    const double &baseline, const double &response, const int &number_trials,
		    const double &variance, double &result){
  result = result - 1/(2*variance) * ( pow(response - level_proposed - baseline, 2.0) -
				       pow(response - level_current  - baseline, 2.0) );
}

/*Likelihood ratio on log scale for Binomial model truncated at 0 with success probability on
  logit scale modelled by monotonic function

INPUT:
 level_proposed - Proposed level of the point
 level_current  - Currently assigned level of the point
 baseline       - Currently assigned baseline level
 response       - Number of successes
 number_trials  - Number of trials
 param          - Empty entry
 result         - Current likelihood difference
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void trBinomial_model(const double &level_proposed, const double &level_current,
		      const double &baseline, const double &response, const int &number_trials,
		      const double &param, double &result){

  double prob_prop = 1 - exp( -number_trials * log(1 + exp(level_proposed + baseline)) ); 
  double prob_curr = 1 - exp( -number_trials * log(1 + exp(level_current  + baseline)) );

  result = result + (level_proposed - level_current) * response - number_trials *
    (log(1+exp(level_proposed + baseline)) - log(1+exp(level_current + baseline))) -
    log(prob_prop) + log(prob_curr);
}

/*Likelihood ratio on log scale for Generalized Pareto model with scale on log scale modelled
  by monotonic function

INPUT:
 level_proposed - Proposed level of the point
 level_current  - Currently assigned level of the point
 baseline       - Currently assigned baseline level
 response       - Observation of response variable
 number_trials  - Empty
 shape          - Shape parameter of the GPD random variable
 result         - Current likelihood difference
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void GPD_model(const double &level_proposed, const double &level_current,
	       const double &baseline, const double &response, const int &number_trials,
	       const double &shape, double &result){

  result = result - level_proposed + level_current - 
    (1/shape + 1) * ( log( 1 + shape * response/exp(level_proposed + baseline) ) -
		      log( 1 + shape * response/exp(level_current  + baseline) ) );
}

/*Likelihood ratio on log scale for Bernoulli model with succes probability on logit scale modelled
  by monotonic function

INPUT:
 level_proposed  - Proposed level of the point
 level_current   - Currently assigned level of the point
 baseline        - Currently assigned baseline level
 response        - 1 if success and 0 otherwise
 number_trials   - Number of trials for 0-1 observations 
 param           - Empty entry
 result          - Current likelihood ratio on log scale
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void Bernoulli_model(const double &level_proposed, const double &level_current,
		     const double &baseline, const double &response,
		     const int &number_trials, const double &param, double &result){
  result = result + 
    response * ( log(1 - exp(-number_trials * log(1 + exp(level_proposed + baseline)))) - 
		 log(1 - exp(-number_trials * log(1 + exp(level_current  + baseline)))) ) +
    (1 - response) * number_trials * ( log(1 + exp(level_current  + baseline)) - 
				       log(1 + exp(level_proposed + baseline)) );
}

/*Likelihood ratio on log scale for Binomial model with success probability on logit scale
 modelled by monotonic function for proposed new baseline

INPUT:
 n                 - Number of observations
 level             - Vector of currently assigned levels to the responses
 baseline_current  - Currently assigned baseline
 baseline_proposed - Proposed baseline
 response          - Observations of the response variable
 number_trials     - Number of trials
 param             - Empty entry
 result            - Current likelihood ratio on log scale
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void Binomial_model_re(const int &n, const vector<double> &level, const double &baseline_current,
		       const double &baseline_proposed, const vector<double> &response,
		       const vector<int> &number_trials, const double &param, double &result){

  for(int i=0; i<n; ++i){
    result = result + (baseline_proposed - baseline_current) * response[i] -
      number_trials[i]*( log(1 + exp(level[i] + baseline_proposed)) -
			 log(1 + exp(level[i] + baseline_current )) );
  }
}

/*Likelihood ratio on log scale for Poisson model with rate on log scale modelled by monotonic
 function for proposed new baseline

INPUT:
 n                 - Number of observations
 level             - Vector of currently assigned levels to the responses
 baseline_current  - Currently assigned baseline
 baseline_proposed - Proposed baseline
 response          - Observations of the response variable
 number_trials     - Number of trials
 param             - Empty entry
 result            - Current likelihood ratio on log scale
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void Poisson_model_re(const int &n, const vector<double> &level, const double &baseline_current,
		      const double &baseline_proposed, const vector<double> &response,
		      const vector<int> &number_trials, const double &param, double &result){

  for(int i=0; i<n; ++i){
    result = result + response[i] * (baseline_proposed - baseline_current) -
      (exp(level[i] + baseline_proposed) -  exp(level[i] + baseline_current));
  }
}

/*Likelihood ratio on log scale for Gaussian model with mean modelled by monotonic function
 for proposed new baseline

INPUT:
 n                 - Number of observations
 level             - Vector of currently assigned levels to the responses
 baseline_current  - Currently assigned baseline
 baseline_proposed - Proposed baseline
 response          - Observation of the response variable
 number_trials     - Empty entry
 variance          - Estimated variance of the Gaussian random variable (constant) 
 result            - Current likelihood ratio on log scale
OUTPUT:
 Updated likelihood ratio of the proposed and current point process
*/

void Gaussian_model_re(const int &n, const vector<double> &level, const double &baseline_current,
		       const double &baseline_proposed, const vector<double> &response,
		       const vector<int> &number_trials, const double &variance, double &result){

  for(int i=0; i<n; ++i){
    result = result - 1/(2*variance) * ( pow(response[i] - level[i] - baseline_proposed, 2) -
					 pow(response[i] - level[i] - baseline_current , 2) );
  }
}

/*Likelihood ratio on log scale for Binomial model truncated at 0 with success probability on
  logit scale modelled by monotonic function for proposed new baseline

INPUT:
 n                 - Number of observations
 level             - Vector of currently assigned levels to the responses
 baseline_current  - Currently assigned baseline
 baseline_proposed - Proposed baseline
 response          - Observations of the response variable
 number_trials     - Number of trials
 param             - Empty entry
 result            - Current likelihood ratio on log scale
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void trBinomial_model_re(const int &n, const vector<double> &level, const double &baseline_current,
			 const double &baseline_proposed, const vector<double> &response,
			 const vector<int> &number_trials, const double &param, double &result){
 
 for(int i=0; i<n; ++i){

   // Derive probabilities for Binomial exceeding 0 for proposed and current levels
   double prob_prop = 1 - exp(-number_trials[i] * log(1 + exp(level[i] + baseline_proposed))); 
   double prob_curr = 1 - exp(-number_trials[i] * log(1 + exp(level[i] + baseline_current )));
   
   //Update likelihood ratio
   result = result - log(prob_prop) + log(prob_curr) + 
     (baseline_proposed - baseline_current) * response[i] - number_trials[i] * 
     ( log(1 + exp(level[i] + baseline_proposed)) - log(1 + exp(level[i] + baseline_current)) );
  }
}

/*Likelihood ratio on log scale for Generalized Pareto model with scale on log scale modelled
  by monotonic function for proposed new baseline

INPUT:
 n                 - Number of observations
 level             - Vector of currently assigned levels to the responses
 baseline_current  - Currently assigned baseline
 baseline_proposed - Proposed baseline
 response          - Observations of response variable
 number_trials     - Empty
 shape             - Shape parameter of the GPD random variable
 result            - Current likelihood ratio on log scale
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void GPD_model_re(const int &n, const vector<double> &level, const double &baseline_current,
		  const double &baseline_proposed, const vector<double> &response,
		  const vector<int> &number_trials, const double &shape, double &result){

  for(int i=0; i<n; ++i){
    result = result - baseline_proposed  + baseline_current - 
      (1/shape + 1) * ( log(1 + shape * response[i]/exp(level[i] + baseline_proposed)) - 
			log(1 + shape * response[i]/exp(level[i] + baseline_current )) );
  }
}

/*Likelihood ratio on log scale for Bernoulli model with succes probability on logit scale modelled
  by monotonic function for proposed new baseline

INPUT:
 n                 - Number of observations 
 level             - Vector of currently assigned levels
 baseline_current  - Currently assigned random effect
 baseline_proposed - Proposed random effect
 response          - Observations of the response variable
 number_trials     - Number of trials
 param             - Empty entry
 result            - Current likelihood ratio on log scale
OUTPUT:
 Updated likelihood ratio of proposed and current point process
*/

void Bernoulli_model_re(const int &n, const vector<double> &level, const double &baseline_current,
			const double &baseline_proposed, const vector<double> &response,
			const vector<int> &number_trials, const double &param, double &result){

  for(int i=0; i<n; ++i){
    result = result + 
      response[i] * (log(1 - exp(-number_trials[i] * log(1 + exp(level[i]+baseline_proposed)))) - 
		     log(1 - exp(-number_trials[i] * log(1 + exp(level[i]+baseline_current)))) ) + 
      (1 - response[i]) * number_trials[i] * ( log(1 + exp(level[i] + baseline_current )) - 
					       log(1 + exp(level[i] + baseline_proposed)) );
  }
}
