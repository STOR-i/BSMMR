#include <vector>
#include <iostream>
#include <stdlib.h>
#include <random>
#include <string>
#include "utilities.hpp"
#include "cross_validation.hpp"
#include "updates.hpp"
#include "analysis.hpp"
#include "modelfunctions.hpp"

using namespace std;

int main(){
  
  //Initialize parameters for the algorithm and read them from the configuration fioe "input.txt"
  string model, set, cross_val;
  int  number_repetitions, number_iterations, burn_in, size_thinning, max_jumps, K, grid_size,
    convergence_samples_size, dim;
  double p_birth, p_death, alpha_param, beta_param, DELTA_MIN, DELTA_MAX, omega, omega_0, p, q,
    eta, tau, alpha_tau, beta_tau;
  read_configurationfile(number_iterations, burn_in, size_thinning, p_birth, p_death, max_jumps,
			 model, alpha_param, beta_param, K, dim, DELTA_MIN, DELTA_MAX, set,
			 cross_val, omega, omega_0, number_repetitions, p, q, eta, tau,
			 alpha_tau, beta_tau, grid_size, convergence_samples_size);

  //Import the data from "data.txt"
  vector<vector<vector<double> > > covariates;
  vector<vector<int> > neighbours, number_trials;
  vector<vector<double> > weights, weights_baseline, assigned_level_curr, response;
  
  import_data(covariates, response, number_trials, neighbours, assigned_level_curr, weights,
	      weights_baseline, dim, DELTA_MIN, model, K);
 
  //Initialize further parameters to perform the algorithm
  vector<double> baseline, sum_weights_baseline, param, mean;
  vector<vector<int> > ind, number_jumps_subprocesses, low_sets, upp_sets;
  vector<vector<double> > sum_levels;
  vector<int> number_jumps;
  vector<vector<vector<int> > > pos_cov;
  vector<vector<vector<vector<double> > > > processes;
  initialise_parameters(processes, ind, baseline, param, low_sets, upp_sets, pos_cov, sum_levels,
			number_jumps_subprocesses, number_jumps, mean, sum_weights_baseline,
			weights_baseline, dim, K); 

  //Derive sample space for each region and derive spaces A_k,k'
  vector<vector<double> > volume_sample_space;
  vector<vector<vector<double> > > sample_boundaries;
  vector<vector<vector<vector<double> > > > A;
  derive_volume(A, sample_boundaries,volume_sample_space,covariates, dim, K, set, neighbours, ind);

  //Initialize the log-likelihood functions for the probability model
  void(*model_log_like)(const double &, const double &, const double &, const double &,
			const int &, const double &, double &);
  void (*model_log_like2)(const int &, const vector<double> &, const double &, const double &,
			  const vector<double> &, const vector<int> &, const double &, double &);

  //Case 1: Binomial model 
  if(model == "Binomial"){ 
    model_log_like  = Binomial_model;
    model_log_like2 = Binomial_model_re;
  }
  //Case 2: Poisson model
  else if(model == "Poisson"){
    model_log_like  = Poisson_model;
    model_log_like2 = Poisson_model_re;
  }
  //Case 3: Gaussian model
  else if(model == "Gaussian"){
    model_log_like  = Gaussian_model;
    model_log_like2 = Gaussian_model_re;
  }
  //Case 4: Binomial model for reponse being strictly positive (trundcated Binomial)
  else if(model == "tBinom"){
    model_log_like  = trBinomial_model;
    model_log_like2 = trBinomial_model_re;
  }
  //Case 5: Generalized Pareto model
  else if(model == "GPD"){
    model_log_like  = GPD_model;
    model_log_like2 = GPD_model_re;
  }
  //Case 6: Bernoulli model for response being 0 or 1
  else if(model == "Bernoulli"){
    model_log_like  = Bernoulli_model;
    model_log_like2 = Bernoulli_model_re;
  }
  //Otherwise, stop algorithm and send message that model is not defined in the current version
  else{
    cout<<"Undefined model"<< endl;
    exit(1);
  }

  //Perform cross validation to find optimal omega - if no omega given
  int number_folds = 10;
  if(cross_val=="true"){
    cout<<"Performing cross validation"<<endl;
    omega = optimal_omega(omega_0, number_repetitions, number_folds, K, model, neighbours, set,
			  weights, sample_boundaries, A, volume_sample_space, DELTA_MIN,
			  DELTA_MAX, p, q, dim, eta, tau, alpha_tau, beta_tau, alpha_param,
			  beta_param, p_birth, p_death, max_jumps, covariates, number_trials,
			  response, weights_baseline, model_log_like, model_log_like2);
  }
  
  //Perform RJMCMC algorithm with chosen or derived omega
  uniform_int_distribution<int> sample_subprocess(0, ind.size() -1);
  cout<<"Start MCMC algorithm"<<endl;
  mt19937 generator(200000);

  for(int i=0; i<number_iterations; ++i){

    //Update the marked point processes for Region 1 to K
    for(int k=0; k<K; ++k){
      int ip = sample_subprocess(generator);
      update_process(k, ip, processes, ind, covariates[k], assigned_level_curr[k], baseline[k], 
		     response[k], number_trials[k], param[k], set, neighbours[k], weights[k],
		     sample_boundaries[k], A[k], low_sets, upp_sets, pos_cov,
		     volume_sample_space[k][ip], DELTA_MIN, DELTA_MAX, omega, dim, p, q, eta,
		     p_birth, p_death, sum_levels[k][ip], number_jumps[k],
		     number_jumps_subprocesses[k][ip], max_jumps, generator, model_log_like);
    }

    //Update baseline levels for Region 1 to K
    //update_baseline(baseline, K, response, number_trials, assigned_level_curr, neighbours,
    //		   weights_baseline, sum_weights_baseline, mean, param, tau, generator,
    //		   model_log_like2);

    //Update the spatial smoothing parameter of the baseline levels
    update_tau(tau, sum_weights_baseline, baseline, weights, neighbours, K, alpha_tau, beta_tau,
	       generator);

    //Update additional model parameter if probability model is of Gaussian or GPD type	
    if(model=="Gaussian" || model=="GPD")
      update_param(model,param, K, assigned_level_curr, baseline, response, alpha_param,
		   beta_param, generator);
    
    //Store current status of the sampled Markov chain after burn-in and with thinning
    if(i % size_thinning == 0 && i > burn_in){
      cout<<i<<endl;
      write_to_file(processes, K, dim, ind, baseline, tau, param, model, number_jumps,
		    sample_boundaries);
    }
  }

  //Analyse the output of the RJMCMC algorithm  
  cout<<"Start analysis of the output"<<endl;
  analysis_and_prediction(model, dim, K, grid_size, DELTA_MIN, convergence_samples_size,
			  sample_boundaries, generator);
  cout<<"Algorithm finished succesfully"<<endl;
}
