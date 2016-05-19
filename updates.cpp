#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>
#include "updates.hpp"
#include "modelfunctions.hpp"
#include "calculate_accept_prob.hpp"
#include "function_evaluation.hpp"

using namespace std;

/*Sample the proposal, calculate the acceptance probability and update process accordingly

INPUT:
 k                       - Indice of the region on which the proposal is made
 ip                      - Indice of the subprocesses proposed to be updated
 processes               - Currently estimated processes for all regions 
 ind                     - Index set describing the considered covariate subsets
 observations            - Covariate observations for region k
 level_current           - Currently assigned levels
 random_effect           - Random effect for region k
 number_claims           - Observations of the response variable for region k 
 number_policies         - Number of trials if Binomial model and empty otherwise
 param                   - Additional parameter if necessary for model(Gaussian and GPD)
 set                     - String indicating whether covariate spaces are varying or not
 neighbours              - Vector containing the indices of the neighbours of region k 
 weights                 - Constants d_k,k' specifying the similarity of region k and neighbours
 sample_boundaries       - Lower and upper corner on which points are proposed
 A_k                     - Left and upper corner of A_k,k' for all neighbours k' of region k
 low_sets                - Indics of the processes on which proposal may put constraint
 upp_sets                - Indics of the processes which may put constraint on proposal
 pos_cov                 - Matrix used to compare points from two different set. Gives the indices 
 volume_sample_space     - Volume of the space on which proposed
 DELTA_MIN               - Smallest possible level
 DELTA_MAX               - Largest possible level   
 omega                   - Spatial smoothing parameter
 dim                     - Dimension of the full covariate set
 p,q                     - Values for the distance measure
 eta                     - Parameter penalizing model complexity
 p_birth, p_death        - Probabilities for a proposed birth and death respectively 
 generator               - Random number generator
 sum_levels              - Sum of the levels in the subprocess on which proposed
 number_jumps            - Total number fo support points describing the monotonic function
 number_jumps_subprocess - Number of jumps in the subprocess
 maximum_jumps           - Highest possible number of points allowed in the subprocess
 model_log_like          - Log-likelihood function of the considered probability model
OUTPUT:
 Updated process and values describing the process           
SUBFUNCTIONS CALLED:
 sample_proposal_birth, accept_prob_birth, sample_proposal_death, accept_prob_death,
 sample_proposal_shift, accept_prob_shift
*/

void update_process(const int &k, const int &ip, 
		    vector<vector<vector<vector<double> > > > &processes,
		    const vector<vector<int> > &ind, const vector<vector<double> > &covariates_k,
		    vector<double> &assigned_level_curr, const double &baseline_k, 
		    const vector<double> &response_k, const vector<int>&number_trials_k,
		    const double &param_k, const string &set, const vector<int> &neighbours,
		    const vector<double> &weights, const vector<vector<double> >&sample_boundaries,
		    const vector<vector<vector<double> > > &A_k,
		    const vector<vector<int> > &low_sets, const vector<vector<int> > &upp_sets,
		    const vector<vector<vector<int> > > &pos_cov,
		    const double &volume_sample_space, const double &DELTA_MIN,
		    const double &DELTA_MAX, const double &omega, const int &dim, const double &p,
		    const double &q, const double &eta, const double &p_birth,
		    const double &p_death, double &sum_levels, int &number_jumps,
		    int &number_jumps_subprocess, const int & maximum_jumps, mt19937 &generator,
		    void (*model_log_like)(const double &, const double &, const double &,
					   const double &, const int &, const double &, double &)){
 
  //Sample whether proposed move is "Birth", "Death" or "Switch"
  string proposed_move;
  uniform_real_distribution<double> standard_uniform(0,1); 
  double temp = standard_uniform(generator);//Sample realization from uniform distribution
  //If the number of points is smaller than maximum number, "Birth" can be proposed
  if(number_jumps_subprocess < maximum_jumps){
    if(temp < p_birth) proposed_move = "Birth";
    //If at least one point is in the process "Death" and "Shift" are possible
    else if (number_jumps_subprocess > 0){
      if(temp < p_birth + p_death) proposed_move = "Death";
      else proposed_move="Shift";
    }
    //Else the process remains unchanged if "Birth" is not proposed  
    else return;
  }
  //If the maximum number of points is reached, "Death" and "Shift are only options
  else if(temp < p_birth) return;//Process remains unchanged
  else if(temp < p_birth  + p_death) proposed_move = "Death";
  else proposed_move="Shift";

  vector<vector<double> > proposal;//Initialize the proposal
  vector<double> assigned_level_prop = assigned_level_curr;//Copy the current levels
 
  //Case 1: Proposed move is "Birth"
  if(proposed_move == "Birth"){
    
    //Sample location and mark of point to add to the process and derive values for proposal ratio
    proposal = sample_proposal_birth(ip, processes[k], ind, sample_boundaries, DELTA_MIN,
				     DELTA_MAX, low_sets[ip], upp_sets[ip], pos_cov, generator);
    proposal.push_back( {volume_sample_space} );
  
    //Calculate acceptance probability
    double result = accept_prob_birth(k, processes, ind, proposal, ip, covariates_k,
				      assigned_level_curr, assigned_level_prop, baseline_k,
				      response_k, number_trials_k, param_k, set, neighbours,
				      weights, sample_boundaries, A_k, low_sets[ip], pos_cov,
				      DELTA_MIN, omega, dim, p, q, eta, p_birth, p_death,
				      model_log_like);

    //Accept or reject proposal and update point process if proposal is accepted
    double u = standard_uniform(generator);
    if(log(u) < result){
      processes[k][ip].insert( processes[k][ip].begin() + (int)proposal[0][0], proposal[2] );
    
      //Update variable values that change with the new point 
      assigned_level_curr = assigned_level_prop;//Update levels assigned to the data points
      sum_levels = sum_levels + proposal[2].back();//Sum of levels for points in the process
      ++number_jumps_subprocess;//Number of jumps in the subprocess ip
      ++number_jumps;//Number of points describing the function for Region k
    }
    else return;
  }

  //Case 2: Proposed move is "Death"
  else if(proposed_move=="Death"){

    //Sample point to be deleted from the process and derive values for the proposal ratio
    proposal = sample_proposal_death(ip, processes[k], ind, DELTA_MIN, DELTA_MAX, low_sets[ip],
				     upp_sets[ip], pos_cov, generator);
    proposal.push_back( {volume_sample_space} );

    //Derive acceptance probability and update level_proposed
    double result = accept_prob_death(k, processes, ind, proposal, ip, covariates_k, 
				      assigned_level_curr, assigned_level_prop,
				      baseline_k, response_k, number_trials_k, param_k, set,
				      neighbours, weights, sample_boundaries, A_k, low_sets[ip],
				      pos_cov, DELTA_MIN, omega, dim, p, q, eta, p_birth, p_death,
				      model_log_like);

    //Accept or reject proposal and update point process if proposal is accepted
    double u = standard_uniform(generator);
    if(log(u) < result){
      //Update variable values that change with the removal of the point
      assigned_level_curr = assigned_level_prop;//Update levels assigned to the data points
      sum_levels = sum_levels - proposal[2].back();//Sum of levels for points in the process
      --number_jumps_subprocess;//Number of points in subprocess ip
      --number_jumps;//Number of of points describing function for Region k
    }
    //If proposal is rejected, the point is attached to the process again
    else processes[k][ip].insert(processes[k][ip].begin() + (int)proposal[0][0], proposal[2]);
  }
  
  //Case 3: Proposed move is "Shift"
  else{
    double temp_sum_levels = sum_levels;//Copy the sum of levels in the subprocess 

    //Sample to be shifted and sample new location and mark for this point   
    proposal = sample_proposal_shift(ip, processes[k], ind, sample_boundaries, DELTA_MIN,
				     DELTA_MAX, low_sets[ip], upp_sets[ip], pos_cov,
				     temp_sum_levels, generator);

    //Calculate acceptance probability and update level_proposed
    double result = accept_prob_shift(k, processes, ind, proposal, ip, covariates_k,
				      assigned_level_curr, assigned_level_prop, baseline_k,
				      response_k, number_trials_k, param_k, set, neighbours,
				      weights, sample_boundaries, A_k, low_sets[ip], pos_cov,
				      DELTA_MIN, omega, dim, p, q, sum_levels, temp_sum_levels,
				      model_log_like);

    //Accept or reject proposal and update point process if proposal is accepted
    double u = standard_uniform(generator);
    if(log(u) < result){
      //Update variable values that change with the removal of the point
      assigned_level_curr = assigned_level_prop;//Updated levels assigned to the data points
      sum_levels = temp_sum_levels;//Sum of the levels in the subprocess
    }
    //Otherwise, replace proposal by the orginal point
    else processes[k][ip][(int)proposal[0][0]] =  proposal[2];
  }
}

/*Update the random effects for all K regions

INPUT:
 baseline        - Currently assigned baselines for the K regions
 K               - Number of considered regions
 response        - Observations for all K regions
 number_trials   - Number of trials for the K regions if the considered model is of Binomial type
 assigned_levels - Currently assigned levels to the response variables
 neighbours      - Neighbourhood structure of the K regions
 weights         - Constants describing the pair-wise similarity of the regions
 sum_weights     - Sum of the weights that are associated to a region
 mean            - Weighted average of the neighbouring regions for each of the K regions
 param           - Additional parameter in the model(Gaussian or GPD)
 tau             - Smoothing parameter for the random effects
 generator       - Random number generator
 model_log_like  - Log-likelihood funtion asscoiated to the model
OUTPUT:
 updated baselines
*/

void update_baseline(vector<double> &baseline, const int&K,const vector<vector<double> > &response,
		     const vector<vector<int> > &number_trials,
		     const vector<vector<double> > &assigned_level,
		     const vector<vector<int> > &neighbours,const vector<vector<double> > &weights,
		     const vector<double> &sum_weights, vector<double> &mean,
		     const vector<double> &param, const double &tau, mt19937 &generator,
		     void (*model_log_like)(const int &, const vector<double> &, const double &,
					    const double &, const vector<double> &,
					    const vector<int> &, const double &, double &)){
 
  //Set sample distribution for updating baselines
  normal_distribution<double> normal(0, 0.05);
  uniform_real_distribution<double> standard_uniform(0,1);

  for(int k=0; k<K; ++k){
    //Sample proposal by Random Walk for the baseline for Region k 
    double proposal = baseline[k] + normal(generator);

    //Derive acceptance probability for proposed baseline
    double accept_prob = 0.0;
    model_log_like(assigned_level[k].size(), assigned_level[k], baseline[k], proposal,
		   response[k], number_trials[k], param[k], accept_prob);
    accept_prob = accept_prob - (tau*sum_weights[k])/2*(pow(proposal - mean[k], 2.0) -
							pow(baseline[k] - mean[k], 2.0));
    
    //Sample whether to update the baseline of Region k
    double u = standard_uniform(generator);
    if(log(u) < accept_prob){
      //Update components of the prior for neighbouring regions
      for(int i=0; i<neighbours[k].size(); ++i){
	int j = neighbours[k][i];
	mean[j] = mean[j] - (baseline[k] - proposal) * weights[k][i] / sum_weights[j];
      }
      baseline[k] = proposal;//Update baseline level
    }
  }
}

/*Update the spatial smoothing parameter for the random effects

INPUT:
 tau            - Smoothing parameter for the random effects
 sum_weights    - Sum of the weights that are associated to a region
 baseline       - Currently assgined baseline to the K regions
 weights        - Constants describing similarity of the regions    
 neighbours     - Neighbourhood strcuture of the regions
 K              - Number of considered regions
 alpha_tau      - Shape parameter of the Gamma prior on eta 
 beta_tau       - Rate parameter of the Gamma prior on eta
 generator      - Random number generator
OUTPUT:
 Updated smoothing parameter tau
*/

void update_tau(double &tau, const vector<double> &sum_weights, const vector<double> &baseline,
		const vector<vector<double> > weights, const vector<vector<int> > neighbours,
		const int &K, const double &alpha_tau, const double &beta_tau, mt19937 &generator){
  
  //Derive parameters of the Gamma posterior for the smoothing parameter tau
  double alpha = K/2.0 + alpha_tau;
  double beta = beta_tau;
  for(int k=0; k<K; ++k){
    beta = beta + pow(baseline[k],2) * sum_weights[k];
    for(int m=0; m<neighbours[k].size(); ++m)
      beta = beta - weights[k][m] * baseline[k] * baseline[neighbours[k][m]];
  }

  //Sample new value for tau from the Gamma posterior
  gamma_distribution<double> sample_eta(alpha, 1.0/beta);
  tau = sample_eta(generator);
}

/*Update additional parameters in the model

INPUT:
 model         - Name of the model considered
 param         - Currently assigned parameter values for the K regions
 K             - Number of regions
 level         - Currently assigned levels to the observations
 random_effect - Currently assigned random effects to the K regions
 number_claims - Observations for the K regions
 alpha_param   - First parameter in the prior distribution
 beta_param    - Second parameter in the prior distribution
OUTPUT:
 Updated parameter values
*/
  
void update_param(const string &model, vector<double> &param, const int &K,
		  const vector<vector<double> > &assigned_level, const vector<double> &baseline,
		  const vector<vector<double> > &response, const double &alpha_param,
		  const double &beta_param, mt19937 &generator){

  //Update variance in a Gaussian setting by sampling from the standard Gamma posterior
  if(model=="Gaussian"){
    for(int k=0; k<K; ++k){
      double alpha = assigned_level[k].size()/2.0 + alpha_param;
      double beta = beta_param;
      for(int i=0; i<assigned_level[k].size(); ++i)
	beta = beta + 0.5 * pow(response[k][i] - assigned_level[k][i] - baseline[k], 2); 
      gamma_distribution<double> sample_param(alpha, 1.0/beta);
      param[k] = 1.0/sample_param(generator);
    }
  }
  //Update shape parameter for Generalized Pareto dsitribution using Random Walk Metropolis
  else{
    normal_distribution<double> normal(0, 0.1);
    uniform_real_distribution<double> standard_uniform(0,1);

    for(int k=0; k<K; ++k){ 
      double proposal = param[k] + normal(generator);
      double accept_prob = 0.0;
      for(int i =0; i<assigned_level[k].size(); ++i){
	accept_prob =
	  -(1/proposal+1)*log( 1+proposal*response[k][i]/(assigned_level[k][i] + baseline[k]) ) + 
	  (1/param[k]+1)* log( 1+param[k]*response[k][i]/(assigned_level[k][i] + baseline[k]) ) -
	  (pow(proposal - alpha_param, 2.0) - pow(param[k] - alpha_param, 2.0))/(2*beta_param);
      }
      double u = standard_uniform(generator);
      if(log(u) < accept_prob) param[k] = proposal;
    }
  }
}

/*Sample proposed support point and associated level for a subprocess

INPUT:
 ip                - Indice of the process on the which the proposal is made
 process           - Currently estimated processes for the region the proposal is made on
 ind               - Index set describing the considered covariate subsets
 sample_boudnaries - Lower an upper corner of the area the function is estimated on
 DELTA_MIN         - Smallest possible level
 DELTA_MAX         - Largest possible level  
 low               - Subspaces on which the proposal may put an upper monotonic constraint
 upp               - Subspaces which may put constraint on the proposal
 pos_cov           - Indices used to compare points defined on different covariate subsets
 generator         - Random number generator

OUTPUT:
 proposal with
 [0] Indice where to insert point if accepted
 [1] Lower and upper bound of the potential level
 [2] Sampled location and mark
*/

vector<vector<double> > sample_proposal_birth(const int &ip,
					      const vector<vector<vector<double> > > &process,
					      const vector<vector<int> > &ind,
					      const vector<vector<double> >&sample_boundaries,
					      const double &DELTA_MIN, const double &DELTA_MAX,
					      const vector<int> &low, const vector<int> &upp,
					      const vector<vector<vector<int > > > &pos_cov,
					      mt19937 &generator){
 
  vector<vector<double> > proposal;//Initialize set of return values  

  //Sample proposed location uniformly over the sample space for point process ip
  vector<double> location; 
  for(int i=0; i<ind[ip].size(); ++i){
    uniform_real_distribution<double> sample_covariate(sample_boundaries[0][ind[ip][i]],
						       sample_boundaries[1][ind[ip][i]]);
    location.push_back(sample_covariate(generator));
  }

  //Find indice at which to insert proposal if accepted and store it
  int pos = 0;
  for(int i = 0; i<process[ip].size(); ++i){
    if(process[ip][pos][0] < location[0]) ++pos;
    else break;
  }
  proposal.push_back( {1.0 * pos} );

  //Derive and store the boundaries of the level 
  proposal.push_back(bounds_level(process, location, ip, DELTA_MIN, DELTA_MAX, low, upp, pos_cov));

  //Sample mark for the location and store the derived point
  uniform_real_distribution<double> sample_mark(proposal[1][0], proposal[1][1]);
  location.push_back(sample_mark(generator));
  proposal.push_back(location);
  return(proposal);
}

/*Sample support point and associated level to be proposed for being removed from the subprocess

INPUT:
 ip                - Indice of the process on the which the proposal is made
 process           - Currently estimated processes for the region the proposal is made on
 ind               - Index set describing the considered covariate subsets
 DELTA_MIN         - Smallest possible level
 DELTA_MAX         - Largest possible level  
 low               - Subspaces on which the proposal may put an upper monotonic constraint
 upp               - Subspaces which may put constraint on the proposal
 pos_cov           - Indices used to compare points defined on different covariate subsets
 generator         - Random number generator
OUTPUT:
 proposal with
 [0] Indice where to insert point again if proposal rejected
 [1] Lower and upper bound of the potential level
 [2] Current point which is removed from the process
*/

vector<vector<double> > sample_proposal_death(const int &ip,
					      vector<vector<vector<double> > > &process,
					      const vector<vector<int> > &ind,
					      const double &DELTA_MIN, const double &DELTA_MAX,
					      const vector<int> &low, const vector<int> &upp,
					      const vector<vector<vector<int > > > &pos_cov, 
					      mt19937 &generator){

  vector<vector<double> > proposal;//Initialize set of return values
 
  //Sample point to be deleted from the process uniformily and store its indice and the point
  uniform_int_distribution<int> sample_point(0, process[ip].size() - 1);
  int ind_point = sample_point(generator);
  proposal.push_back( {1.0*ind_point} );
  proposal.push_back( {process[ip][ind_point]} );

  //Remove sampled from the process and derive upper and lower bounds at that point
  process[ip].erase(process[ip].begin() + ind_point);
  proposal.insert( proposal.begin() + 1, bounds_level(process, proposal[1], ip, DELTA_MIN,
						      DELTA_MAX, low, upp, pos_cov) );
  return(proposal);
}

/*Sample proposed support point and associated level for a subprocess for a shift being proposed

INPUT:
 ip                - Indice of the process on the which the proposal is made
 process           - Currently estimated processes for the region the proposal is made on
 ind               - Index set describing the considered covariate subsets
 DELTA_MIN         - Smallest possible level
 DELTA_MAX         - Largest possible level  
 low               - Subspaces on which the proposal may put an upper monotonic constraint
 upp               - Subspaces which may put constraint on the proposal
 pos_cov           - Indices used to compare points defined on different covariate subsets
 sum_levels        - Sum of the levels of the points in the subprocess
 generator         - Random number generator

OUTPUT:
 proposal with
 [0] Indice of the point to be shifted
 [1] Lower and upper bound of the potential level
 [2] Current point which is removed from the process
 [3] Proposed support point and level 
*/

vector<vector<double> > sample_proposal_shift(const int &ip,
					      vector<vector<vector<double> > > &process,
					      const vector<vector<int> > &ind,
					      const vector<vector<double> > &sample_boundaries,
					      const double &DELTA_MIN, const double &DELTA_MAX,
					      const vector<int> &low, const vector<int> &upp,
					      const vector<vector<vector<int > > > &pos_cov,
					      double &sum_levels, mt19937 &generator){

  vector<vector<double> > proposal;
  //Sample point to be shifted with probability proportional to its level
  uniform_real_distribution<double> standard_uniform(0,1);
  double u = standard_uniform(generator);
  double u_temp = (process[ip][0].back() - DELTA_MIN)/sum_levels;
  int idp = 0;//Initialize index of the point to be shifted 
  while(u_temp  < u){
    ++idp;
    u_temp = u_temp + (process[ip][idp].back()-DELTA_MIN)/sum_levels;  
  }

  //Store indice of the point to be shifted and current location and mark
  proposal.push_back( {1.0*idp} );
  proposal.push_back( process[ip][idp] );
  process[ip].erase(process[ip].begin() + idp);//Remove point from the process
  //Derive lower and upper boundaries and sample new location 
  vector<vector<double> > boundaries = bounds_location(ip, proposal[1], process, ind,
						       sample_boundaries, low, upp, pos_cov);
  vector<double> location;
  for(int i=0; i<boundaries.size(); ++i){
    uniform_real_distribution<double> sample_covariate(boundaries[i][0], boundaries[i][1]);
    location.push_back(sample_covariate(generator));
  }
  //Derive lower and upper boundaries for the level at the proposed location and store them
  vector<double> boundaries_level = bounds_level(process, location, ip,  DELTA_MIN,
						 DELTA_MAX, low, upp, pos_cov);
  proposal.insert(proposal.begin() + 1, boundaries_level);
  //Sample mark and store proposal as well as update sum of levels
  uniform_real_distribution<double> sample_mark(boundaries_level[0], boundaries_level[1]);
  location.push_back(sample_mark(generator));
  proposal.push_back(location);
  sum_levels = sum_levels - proposal[2].back() + proposal[3].back();
  return(proposal);
}
 
