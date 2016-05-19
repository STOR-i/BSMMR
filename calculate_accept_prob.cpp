#include <vector>
#include <iostream>
#include <cmath> 
#include <stdlib.h>
#include <iomanip>
#include <algorithm>
#include "calculate_accept_prob.hpp"
#include "function_evaluation.hpp"
#include "calculate_integral_birth.hpp"
#include "modelfunctions.hpp"

using namespace std;

/*Calculation of the acceptance probability for a proposed "Birth"

INPUT:
 k                   - Indice of the region on which the proposal is made
 processes           - Currently estimated processes for all regions 
 ind                 - Index set describing the considered covariate subsets
 ip                  - Indice of the point process on which "Birth" is proposed 
 covariates_k        - Covariate observations of the data for Region k
 assigned_level_curr - Currently assigned levels
 assigned_level_prop - Proposed assigned levels to be updated
 baseline_k          - Assigned baseline level for Region k
 response_k          - Observations of the response variable for Region k
 number_trials_k     - Number of trials if Binomial model and empty otherwise forRegion k
 param_k             - Additional parameter if necessary for model (Gaussian and GPD)
 set                 - String indicating whether covariate spaces are varying or not
 neighbours          - Vector containing the indices of the neighbours of region k 
 weights             - Constants d_k,k' specifying the similarity of Region k and its neighbours
 sample_boundaries   - Lower and upper corner on which points are proposed
 A_k                 - Left and upper corner of A_k,k' for all neighbours k' of Region k
 low_sets            - Containing indices of the processes on which proposal may put constraint 
 pos_cov             - Indices to be compared if points are defined on different subsets
 DELTA_MIN           - Smallest possible level
 omega               - Spatial smoothing parameter
 dim                 - Dimension of the full covariate set
 p,q                 - Values for the distance measure
 model_log_like      - Log-likelihood function of the considered probability model
OUTPUT:
 result - acceptance probability for the proposed birth            
SUBFUNCTIONS CALLED:
 "calculate_integral_birth.cpp":
  integral_birth 
 "function_evaluation.cpp":
  check_monotonicity 
 "modelfunctions.cpp":
  model_log_like
*/

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
						const double &, double &)){
			 
  double result = 0.0;

  //STEP 1: Calculate likelihood ratio and update assigned levels for the proposal
  for(int i=0; i<assigned_level_curr.size(); ++i){

    //CASE 1: If current level higher than proposed one, assigned level won't change
    if(proposal[2].back() <= assigned_level_curr[i]) continue;

    //CASE 2: If point puts constraint on proposal, update level and likelihood ratio
    else if( check_monotonicity(proposal[2], covariates_k[i], ind[ip], 0) == true){
      assigned_level_prop[i] = proposal[2].back();
      model_log_like(assigned_level_prop[i], assigned_level_curr[i], baseline_k,
		     response_k[i], number_trials_k[i], param_k, result);
      }	
    else continue;
  }

  //STEP 2: Derive ratio for the spatial smoothing component of the prior densities 
  result = result  - omega * integral_birth(set, k, proposal[2], ip, processes, ind, neighbours,
					    weights, dim, DELTA_MIN, sample_boundaries, A_k,
					    low_sets, pos_cov, proposal[1][0], p, q);

  //STEP 3: Calculate proposal ratio and ratio in model complexity
  result = result  + log(proposal[1][1] - proposal[1][0]) + log(proposal[3][0]) - 
    log(processes[k][ip].size() + 1)  + log(p_death) - log(p_birth) + log(1-1/eta);
  return(result);
}

/*Calculation of the acceptance probability for a proposed "Death"

INPUT:
 k                   - Indice of the region on which the proposal is made
 processes           - Currently estimated processes for all regions 
 ind                 - Index set describing the considered covariate subsets
 ip                  - Indice of the point process on which "Death" is proposed 
 covariates_k        - Covariate observations of the data for Region k
 assigned_level_curr - Currently assigned levels
 assigned_level_prop - Proposed assigned levels to be updated
 baseline_k          - Assigned baseline level for Region k
 response_k          - Observations of the response variable for Region k
 number_trials_k     - Number of trials if Binomial model and empty otherwise forRegion k
 param_k             - Additional parameter if necessary for model (Gaussian and GPD)
 set                 - String indicating whether covariate spaces are varying or not
 neighbours          - Vector containing the indices of the neighbours of region k 
 weights             - Constants d_k,k' specifying the similarity of Region k and its neighbours
 sample_boundaries   - Lower and upper corner on which points are proposed
 A_k                 - Left and upper corner of A_k,k' for all neighbours k' of Region k
 low_sets            - Containing indices of the processes on which proposal may put constraint 
 pos_cov             - Indices to be compared if points are defined on different subsets
 DELTA_MIN           - Smallest possible level
 omega               - Spatial smoothing parameter
 dim                 - Dimension of the full covariate set
 p,q                 - Values for the distance measure
 model_log_like      - Log-likelihood function of the considered probability model
OUTPUT:
 result - Acceptance probability for the proposed death           
SUBFUNCTIONS CALLED:
 "calculate_integral_birth.cpp":
  integral_birth 
 "function_evaluation.cpp":
  function_evaluation
 "modelfunctions.cpp":
  model_log_like
*/

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
						const int &, const double &, double &)){

  double result = 0.0;

  //STEP 1: Calculate likelihood ratio and update assigned levels for the proposal
  for(int i=0; i<assigned_level_curr.size(); ++i){

    //If assigned level is equal to the one to be removed, derive new level and likelihood ratio
    if( abs(assigned_level_curr[i]-proposal[2].back()) < 1e-10){       
      assigned_level_prop[i] = function_evaluation(processes[k], ind, dim, covariates_k[i],
						   proposal[1][0]);
      model_log_like(assigned_level_prop[i], assigned_level_curr[i], baseline_k,
		     response_k[i], number_trials_k[i], param_k, result);
    }
    else continue;//Otherwise, level won't change
  }
  
  //STEP 2: Derive ratio of the prior densities on log scale
  result = result  + omega * integral_birth(set, k, proposal[2], ip, processes, ind, neighbours,
					    weights, dim, DELTA_MIN, sample_boundaries, A_k,
					    low_sets, pos_cov, proposal[1][0], p, q);

  //STEP 3: Calculate ratio in proposals and model complexity
  result = result  - log(proposal[1][1] - proposal[1][0]) - log(proposal[3][0]) + 
    log(processes[k][ip].size()+1) - log(p_death) + log(p_birth) - log(1-1/eta);
  return(result);
}

/*Calculation of the acceptance probability for a proposed "Shift"

INPUT:
 k                   - Indice of the region on which the proposal is made
 processes           - Currently estimated processes for all regions 
 ind                 - Index set describing the considered covariate subsets
 ip                  - Indice of the point process on which "Shift" is proposed 
 covariates_k        - Covariate observations of the data for Region k
 assigned_level_curr - Currently assigned levels
 assigned_level_prop - Proposed assigned levels to be updated
 baseline_k          - Assigned baseline level for Region k
 response_k          - Observations of the response variable for Region k
 number_trials_k     - Number of trials if Binomial model and empty otherwise forRegion k
 param_k             - Additional parameter if necessary for model (Gaussian and GPD)
 set                 - String indicating whether covariate spaces are varying or not
 neighbours          - Vector containing the indices of the neighbours of region k 
 weights             - Constants d_k,k' specifying the similarity of Region k and its neighbours
 sample_boundaries   - Lower and upper corner on which points are proposed
 A_k                 - Left and upper corner of A_k,k' for all neighbours k' of Region k
 low_sets            - Containing indices of the processes on which proposal may put constraint 
 pos_cov             - Indices to be compared if points are defined on different subsets
 DELTA_MIN           - Smallest possible level
 omega               - Spatial smoothing parameter
 dim                 - Dimension of the full covariate set
 p,q                 - Values for the distance measure
 sum_levels          - Sum of the levels in the current subprocess    
 temp_sum_levels     - Sum of the levels in the subprocess if the proposal would be accepted
 model_log_like      - Log-likelihood function of the considered probability model
OUTPUT:
 result - Acceptance probability for the proposed birth switch
SUBFUNCTIONS CALLED:
  calculate_likelihood_shift
 "calculate_integral_birth.cpp":
  integral_birth 
 "modelfunctions.cpp":
  model_log_like
*/

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
			 const double &q, const double &sum_levels, const double &temp_sum_levels,
			 void (*model_log_like)(const double &, const double &,
						const double &,	const double &, const int&,
						const double &, double &)){
			  
  double result = 0.0;
 
  //STEP 1: Derive ratio in the prior densities
  //STEP 1(a): Treat current point as to be removed from the process, i.e. calculate like "Death"
  result = result + omega * integral_birth(set, k, proposal[2], ip, processes, ind, neighbours,
					   weights, dim, DELTA_MIN, sample_boundaries, A_k,
					   low_sets, pos_cov, proposal[1][0], p, q);

  //STEP 1(b): Treat proposed point as to be removed from the process, i.e. calculate like "Birth"
  result = result - omega * integral_birth(set, k, proposal[3], ip, processes, ind, neighbours,
					   weights, dim, DELTA_MIN, sample_boundaries, A_k,
					   low_sets, pos_cov, proposal[1][0], p, q);  

  //STEP 2: Calculate likelihood ratio after first attaching proposed point to the process
  processes[k][ip].insert(processes[k][ip].begin() + (int)proposal[0][0], proposal[3]);
  calculate_likelihood_shift(result, proposal, ip, processes[k], ind, dim, baseline_k,
			     covariates_k, assigned_level_curr, assigned_level_prop, 
			     response_k, number_trials_k, param_k, model_log_like);

  //Add the quotient of the proposal density to the acceptance probability
  result = result + log(proposal[3].back() - DELTA_MIN) - log(proposal[2].back() - DELTA_MIN) -
    log(temp_sum_levels) + log(sum_levels);
  return result;
}  

/*Derive likelihood ratio on log-scale for a proposed "Shift" and update levels for the data points

INPUT: 
 result              - Return variable to which the likelihood ratio is added
 proposal            - Information describing the proposal
 ip                  - Indice of the point process on which "Shift" is proposed 
 process             - Point processes describing the proposed monotonic function
 ind                 - Indices of the covariate subsets on which processes are defined 
 covariates_k        - Covariate observations of the data for Region k
 assigned_level_curr - Currently assigned levels
 assigned_level_prop - Proposed assigned levels to be updated
 baseline_k          - Assigned baseline level for Region k
 response_k          - Observations of the response variable for Region k
 number_trials_k     - Number of trials if Binomial model and empty otherwise forRegion k
 param_k             - Additional parameter if necessary for model (Gaussian and GPD)
 model_log_like      - Log-Likelihood function for the considered model
OUTPUT:
 Derived ratio of the likelihood functions and updated assigned levels for the data points 
SUBFUNCTIONS CALLED:
 "function_evaluation.cpp":
  check_monotonicity, function_evaluation 
 "modelfunctions.cpp":
  model_log_like
*/

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
						      const double &, double &)){
  
  //Case 1: The proposed level is higher than the current level
  if(proposal[2].back() < proposal[3].back()){
    
    for(int i=0; i<assigned_level_curr.size(); ++i){
      
      //Case a: If currently assigned level is higher than proposal, assigned level won't change 
      if(assigned_level_curr[i] > proposal[3].back()) continue;
      
      //Case b: If proposal puts constraint on point, update assigned level and likelihood ratio  
      else if(check_monotonicity(proposal[3], covariates_k[i], ind[ip], 0) == true){
	assigned_level_prop[i] = proposal[3].back();
	model_log_like(assigned_level_prop[i], assigned_level_curr[i], baseline_k,
		       response_k[i], number_trials_k[i], param_k, result);
      }
      
      //Case c: If currnet level equal to the point to be shifted, derive level and update ratio 
      else if( abs(assigned_level_curr[i] - proposal[2].back()) < 1e-10){
	assigned_level_prop[i] = function_evaluation(process, ind, dim, covariates_k[i],
						     proposal[1][0]);
	model_log_like(assigned_level_prop[i], assigned_level_curr[i], baseline_k,
		       response_k[i], number_trials_k[i], param_k, result);
      }
      else continue;
    }
  }
  
  //Case 2: The proposed level is smaller than current level
  else{
    
    for(int i=0; i<assigned_level_curr.size(); ++i){
      
      //Case a: If currnet level equal to the point to be shifted, derive level and update ratio
      if( abs(proposal[2].back()-assigned_level_curr[i]) < 1e-10){
	assigned_level_prop[i] = function_evaluation(process, ind, dim, covariates_k[i],
						     proposal[1][0]);
	model_log_like(assigned_level_prop[i], assigned_level_curr[i], baseline_k,
		       response_k[i], number_trials_k[i], param_k, result);
      }
      
      //Case b: If currently assigned level is higher than proposal, assigned level won't change 
      else if(assigned_level_curr[i] > proposal[3].back()) continue;
      
      //Case c: If proposal puts constraint on point, update assigned level and likelihood ratio 
      else if (check_monotonicity(proposal[3], covariates_k[i], ind[ip], 0) == true){
	assigned_level_prop[i] = proposal[3].back();
	model_log_like(assigned_level_prop[i], assigned_level_curr[i], baseline_k,
		       response_k[i], number_trials_k[i], param_k, result);
      }
      else continue;
    }
  }
}

/*Calculation of the acceptance probability for a proposed baseline

INPUT:
 assigned_level_curr - Currently assigned levels to the observations
 reponse_k           - Observations of the response variable
 baseline_proposed   - Proposed random effect for the region
 baseline_current    - Current random effect for the region
 number_trials_k     - Number of trials(if Binomial model) and empty otherwise
 param_k             - Additional parameter if necessary for model (Gaussian and GPD)
 mean                - Weighted average of the random effects of the neighbouring regions
 sum_weights         - Sum of the weights associated to the region on which the proposal is made
 rho                 - Spatial moothing parameter for the random effects
 model_log_like      - Log-Likelihood function for the considered model
OUTPUT:
 acceptance probability for the proposed baseline
SUBFUNCTIONS CALLED:
 "modelfunctions.cpp":
  model_log_like
*/

double accept_prob_baseline(const vector<double> &assigned_level_curr,
			    const vector<double> &response_k, const double &baseline_prop,
			    const double &baseline_curr, const vector<int> &number_trials_k,
			    const double &param_k, const double &mean, const double &sum_weights,
			    const double &rho,
			    void (*model_log_like)(const int &, const vector<double> &,
						   const double &, const double &,
						   const vector<double> &, const vector<int> &,
						   const double &, double &)){
  
  double result = 0.0;

  //Update likelihood ratio
  model_log_like(assigned_level_curr.size(), assigned_level_curr, baseline_curr,  baseline_prop,
		 response_k, number_trials_k, param_k, result);

  //Derive quotient of the prior densities
  result = result - (rho*sum_weights)/2*(pow(baseline_prop - mean, 2.0) -
					 pow(baseline_curr - mean, 2.0));
  return(result);
}
