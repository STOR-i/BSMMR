#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <fstream>
#include <stdlib.h>
#include "cross_validation.hpp"
#include "utilities.hpp"
#include "updates.hpp"

using namespace std;

/*Derive appropriate value for omega by efficient global optimization and m-fold cross validation

INPUT:
 omega_0             - Inital value for omega 
 rep                 - Number of repetitions of the m-fold cross validation
 number_folds        - Number of folds in which the data set is split
 K                   - Number of regions
 model               - Probability model
 neighbours          - Neighbourhood structure of the regions
 set                 - Indicating whether covariate spaces are varying or not
 weights             - Constants d_k,k' specifying the similarity of region k and its neighbours
 sample_boundaries   - Lower left and upper right corner of the sample space for proposals
 A                   - Left and upper corner of A_k,k' for all apirs k and k'
 volume_sample_space - Volume of the space on which proposed
 DELTA_MIN           - Smallest possible level
 DELTA_MAX           - Largest possible level  
 dim                 - Dimension of the full covariate space
 p,q                 - Values for the distance measure
 eta                 - Parameter penalizing model complexity
 alpha_tau, beta_tau - Hyperparameters for the prior for spatial smoothing of baselines 
 p_birth, p_death    - Probabilities for a proposed birth and death respectively 
 max_jumps           - Highest possible number of points allowed in a point process
 covariates          - Covariate observations for all regions
 response            - Observations of the response variable for the regions 
 number_trials       - Number of trials if Binomial model and empty otherwise
 weights_baseline    - Constants describing the similarity of the regions for the baseline
 model_log_like      - Log-likelihood function for one observation
 model_log_like2     - Log-likelihood function for several observations
OUTPUT:
 optimal_omega - Omega value giving the smallest squared error of observed and predicted output 
SUBFUNCTIONS CALLED:
 cross_validation, BayesOptim.R
*/

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
 					     const vector<int> &, const double &, double &)){
					
  //Set vector of considered values for omega with two entries and initialize remaining variables  
  vector<double> omega( {0, omega_0} ); 
  vector<double> MSE; 
  double minimum_MSE;
  double optimal_omega = 0.0;
  ifstream proposed_omega;//Set filestream from which to read next proposal for omega

  //Perform cross validation for omega = 0 and store result as smallest MSE found so far 
  MSE.push_back(cross_validation(omega[0], rep, number_folds, K, model, neighbours, set,  weights,
				 sample_boundaries, A, volume_sample_space, DELTA_MIN, DELTA_MAX,
				 p, q, dim, eta, tau, alpha_tau, beta_tau, alpha_param, beta_param,
				 p_birth, p_death, max_jumps, covariates, number_trials,
				 response, weights_baseline, model_log_like, model_log_like2)); 
  minimum_MSE = MSE[0];
  
  //Increase initial value for omega until MSE is much larger than MSE for no spatial smoothing
  while(MSE.back() < 1.10*MSE[0]){
    //Perform cross validation for the latest omega
    cout<<omega.back()<<endl;
    MSE.push_back(cross_validation(omega.back(), rep, number_folds, K, model, neighbours, set,
				   weights, sample_boundaries, A, volume_sample_space, DELTA_MIN,
				   DELTA_MAX, p, q, dim, eta, tau, alpha_tau, beta_tau,
				   alpha_param, beta_param, p_birth, p_death, max_jumps,
				   covariates, number_trials, response, weights_baseline,
				   model_log_like, model_log_like2));
 
    //If derived MSE fullfills criterion set a new value for omega to be considered
    if(MSE.back() > 1.10 * MSE[0])
      omega.push_back( 0.5*omega.back() );
    //Otherwise increase last considered value for omega by factor 10 and update detected optimium
    else{
      if(MSE.back() < minimum_MSE){
	minimum_MSE = MSE.back();
	optimal_omega = omega.back();
      }
      omega.push_back(10 * omega.back());
    }
  }

  /*Apply Bayesian optimization using the derived values and upper bound until the expected
    improvement falls below critical value*/
  double critical_value = 0.0001;
  double omega_new, expected_improvement; 

  while(minimum_MSE!=0){//Loop without condition to stop
    cout<<omega.back()<<endl;
    //Perform cross validation for latest proposal for omega
    MSE.push_back(cross_validation(omega.back(), rep, number_folds, K, model, neighbours, set,
				   weights, sample_boundaries, A, volume_sample_space, DELTA_MIN,
				   DELTA_MAX, p, q, dim, eta, tau, alpha_tau, beta_tau, 
				   alpha_param, beta_param, p_birth, p_death, max_jumps,
				   covariates, number_trials, response, weights_baseline,
				   model_log_like, model_log_like2));
  
    //Update optimum solution if necessary
    if(MSE.back() < minimum_MSE){
      minimum_MSE = MSE.back();
      optimal_omega = omega.back();
    }
    //Derive proposal for omega and the associated expected improvement using external R file
    int temp = system("Rscript BayesOptim.R");

    //Read proposal for omega and associated expected improvement
    proposed_omega.open("./Results/proposal.txt");
    proposed_omega >> omega_new;
    proposed_omega >> expected_improvement;
    proposed_omega.close();

    //Check whether expected improvement falls below critical level and stop if so
    if(expected_improvement < critical_value * minimum_MSE)
      return optimal_omega;
    //Otherwise, store the proposal for omega and perform next cross-validation
    else omega.push_back(omega_new);
  }
}

/*Perform m-fold cross validation for a given omega and a number of repetitions

INPUT:
 omega               - Proposal for the spatial smoothing for which cross validation is performed
 rep                 - Number of repetitions of the m-fold cross validation
 number_folds        - Number of folds in which the data set is split
 K                   - Number of regions
 model               - Probability model
 neighbours          - Neighbourhood structure of the regions
 set                 - Indicating whether covariate spaces are varying or not
 weights             - Constants d_k,k' specifying the similarity of region k and its neighbours
 sample_boundaries   - Lower left and upper right corner of the sample space for proposals
 A                   - Left and upper corner of A_k,k' for all apirs k and k'
 volume_sample_space - Volume of the space on which proposed
 DELTA_MIN           - Smallest possible level
 DELTA_MAX           - Largest possible level  
 dim                 - Dimension of the full covariate space
 p,q                 - Values for the distance measure
 eta                 - Parameter penalizing model complexity
 alpha_tau, beta_tau - Hyperparameters for the prior for spatial smoothing of baselines 
 p_birth, p_death    - Probabilities for a proposed birth and death respectively 
 max_jumps           - Highest possible number of points allowed in a point process
 covariates          - Covariate observations for all regions
 response            - Observations of the response variable for the regions 
 number_trials       - Number of trials if Binomial model and empty otherwise
 weights_baseline    - Constants describing the similarity of the regions for the baseline
 model_log_like      - Log-likelihood function for one observation
 model_log_like2     - Log-likelihood function for several observations
OUTPUT:
 MSE - mean squared error of predictive mean and observation for all repetitons, folds and regions
SUBFUNCTIONS CALLED:
  generate_cv_data, derive_expectation
 "utilities.cpp"
  initalize_parameters, , 
 "updates.cpp"
  update_process, update_baseline, update_param, update_tau
*/

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
						const vector<int> &, const double &, double &)){

  vector<double> MSE(rep, 0.0);	
  ofstream cross_validation_file;	

  for(int run=0; run<rep; ++run){
    cout<<"Run "<<run+1<<" of "<<rep<<endl;
    for(int fold=0; fold<number_folds; ++fold){
      //Initialize parameters to run the RJMCMC algorithm for the considered omega
      vector<vector<vector<vector<double> > > > processes;
      vector<vector<int> > ind, low_sets, upp_sets, number_jumps_subprocesses; 
      vector<double> baseline, param, mean, sum_weights_baseline; 
      vector<int> number_jumps;
      vector<vector<vector<int> > > pos_cov;
      vector<vector<double> > sum_levels; 
      initialise_parameters(processes, ind, baseline, param, low_sets, upp_sets, pos_cov,
			    sum_levels, number_jumps_subprocesses, number_jumps, mean,
			    sum_weights_baseline, weights_baseline, dim, K);

      //Split the data into test and training data to perform cross validation
      vector<vector<vector<double> > > test_covariates, training_covariates;
      vector<vector<double> > test_response, training_response; 
      vector<vector<int> > test_number_trials, training_number_trials;
      generate_cv_data(test_covariates, test_number_trials, test_response, training_covariates,
		       training_number_trials, training_response, covariates, response,
		       number_trials, run, fold, K, number_folds, model);
      vector<vector<double> > assigned_level_curr;
      vector<vector<vector<double> > > prediction(K, vector<vector<double> > () );
      for(int k=0; k<K; ++k){
	assigned_level_curr.push_back(vector<double> (training_covariates[k].size(), DELTA_MIN));
	prediction[k] = vector<vector<double> > (test_covariates[k].size(), vector<double> () );
      }
      double temp_tau = tau;

      int size_test_data = test_response[0].size();
      for(int k=1; k<K; ++k) size_test_data = test_response[k].size() + size_test_data;	

      //Run RJMCMC algorithm for a sufficient number of iteration steps
      mt19937 generator(100 * run + 10 * fold);
      uniform_int_distribution<int> sample_subprocess(0, ind.size()-1);

      for(int i=0; i<75001; ++i){
	for(int k=0; k<K; ++k){
	  int ip = sample_subprocess(generator);
	  //Update the marked point process ip for region k
	  update_process(k, ip, processes, ind, training_covariates[k], assigned_level_curr[k],
			 baseline[k], training_response[k], training_number_trials[k], param[k],
			 set, neighbours[k], weights[k], sample_boundaries[k], A[k], low_sets,
			 upp_sets, pos_cov, volume_sample_space[k][ip], DELTA_MIN, DELTA_MAX,
			 omega, dim, p, q, eta, p_birth, p_death, sum_levels[k][ip],
			 number_jumps[k], number_jumps_subprocesses[k][ip], max_jumps,
			 generator, model_log_like);
	}
	
	//Update the baselines for the K regions
	//update_baseline(baseline, K, training_response, training_number_trials,
	//		assigned_level_curr, neighbours, weights_baseline, sum_weights_baseline,
	//		mean, param, temp_tau, generator, model_log_like2);

	//Update the spatial smoothing parameter of the random effects
	update_tau(temp_tau, sum_weights_baseline, baseline, weights, neighbours, K,
		   alpha_tau, beta_tau, generator);
	//Update addiotnal parameters if model is Gaussian or Pareto distributed
	if(model=="Gaussian" || model=="GPD")
	  update_param(model,param, K, assigned_level_curr, baseline, training_response,
		       alpha_param, beta_param, generator);

	//Store currently assigned level for the test data after appropriate burn-in period  
	if(i % 100 == 0 && i > 25000){
	  for(int k=0; k<K; ++k){
	    vector<double> current_level =
	      function_evaluation_set(processes[k], ind, dim, test_covariates[k], DELTA_MIN);
	    //Store the current level for each point separately
	    for(int j=0; j<prediction[k].size(); ++j)
	      prediction[k][j].push_back(current_level[j] + baseline[k]);
	  }
	}  
      }

      //Derive mean squared error of the posterior distribution and the true values
      for(int k=0; k<K; ++k){
	vector<double> expectation;
	derive_expectation(model, prediction[k], param[k], test_number_trials[k], expectation);
	for(int j=0; j<expectation.size(); ++j)
	  MSE[run] = MSE[run] + pow(test_response[k][j] - expectation[j], 2)/size_test_data;
      }
    }
  } 

  //Store omega and the derived MSE for each repetition in external .txt file
  double result = 0.0;	
  cross_validation_file.open("./Results/cross_validation.txt", ios::app);
  cross_validation_file<<omega<<" ";
  for(int i=0; i<rep; ++i){
   cross_validation_file<< MSE[i]/number_folds<<" ";
   result = result + MSE[i]/number_folds;
  }
  cross_validation_file <<"\n";
  cross_validation_file.close();

  return result/rep;//Return MSE over all repetitions
} 

/*Split of the original data into training and test data for cross-validation

INPUT:
 response      - Observations fo the repsonse variable for the K regions
 number_trials - Number of trials if the considered model is of Binomial type and empty otherwise
 covariates    - Covariate observations for the K regions
 run           - Indice of the repetition of the cross validation  
 fold          - Number of the fold considered in the run
 K             - Number of regions
 number_folds  - Number of splits considered in each repetition
 model         - Considered probability model
OUTPUT:
 Training and test data for the fold of the cross-validation stored in
 test_covariates, test_number_trials, test_response, training_covariates, training_number_trials
 and training_response
*/

void generate_cv_data(vector<vector<vector<double> > > &test_covariates,
		      vector<vector<int> > &test_number_trials,
		      vector<vector<double> > &test_response,
		      vector<vector<vector<double> > > &training_covariates,
		      vector<vector<int> > & training_number_trials,
		      vector<vector<double> > &training_response,
		      const vector<vector<vector<double> > > &covariates,
		      const vector<vector<double> > &response,
		      const vector<vector<int> > &number_trials, const int &run, const int &fold,
		      const int &K, const int &number_folds, const string &model){

  mt19937 generator(200 * ( run + 1 ) + 1);//Set random number generator in dependency on run

  //Sample sequence of indices for each region separately 
  vector<vector<int> > indices;
  for(int k=0; k<K; ++k){
    vector<int> temp_indices (covariates[k].size(), 0);
    iota(temp_indices.begin(), temp_indices.end(), 0);
    shuffle(temp_indices.begin(), temp_indices.end(), generator);
    indices.push_back(temp_indices);
  }

  //Initliaze subsets of the test and training data
  test_covariates = vector<vector<vector<double> > > (K, vector<vector<double> > () ); 
  test_number_trials = vector<vector<int> > (K, vector<int> () ); 
  test_response = vector<vector<double> > (K, vector<double> () ); 
  training_covariates = vector<vector<vector<double> > > (K, vector<vector<double> > () ); 
  training_number_trials = vector<vector<int> > (K, vector<int> () ); 
  training_response = vector<vector<double> > (K, vector<double> () ); 

  //Split data based on the shuffled indices generated previously
  
  //If model of Binomial type, number_trials has to be split up too
  if(model =="Binomial" || model=="tBinom"){    
    for(int k=0; k<K; ++k){
      //Define the sequence of indices belonging to that fold
      vector<int> temp_indices;
      int temp1 = (int)( 1.0/number_folds * fold * covariates[k].size() );
      int temp2 = (int)( 1.0/number_folds * (fold+1) * covariates[k].size() );
      for(int i=temp1; i<temp2; ++i) temp_indices.push_back(indices[k][i]);
    
      //Loop over the indices as assign data points either to the test or the training data set
      for(int i=0; i<covariates[k].size(); ++i){
	if(find(temp_indices.begin(),temp_indices.end(), i)!=temp_indices.end()){
	  test_number_trials[k].push_back(number_trials[k][i]);
	  test_response[k].push_back(response[k][i]);
	  test_covariates[k].push_back(covariates[k][i]);
	}
	else{
	  training_number_trials[k].push_back(number_trials[k][i]);
	  training_response[k].push_back(response[k][i]);
	  training_covariates[k].push_back(covariates[k][i]);
	}
      }
    }
  }
  else{
    for(int k=0; k<K; ++k){	
      //Define the sequence of indices belonging to that fold
      vector<int> temp_indices;
      int temp1 = (int)( 1.0/number_folds * fold * covariates[k].size() );
      int temp2 = (int)( 1.0/number_folds * (fold+1) * covariates[k].size() );
      for(int i=temp1; i<temp2; ++i) temp_indices.push_back(indices[k][i]);
  
      //Loop over the indices as assign data points either to the test or the training data set
      for(int i=0; i<covariates[k].size(); ++i){
	if(find(temp_indices.begin(),temp_indices.end(), i)!=temp_indices.end()){	  
	  test_response[k].push_back(response[k][i]);
	  test_covariates[k].push_back(covariates[k][i]);
	}
	else{	  
	  training_response[k].push_back(response[k][i]);
	  training_covariates[k].push_back(covariates[k][i]);
	}
      }
    } 
  }
}

/*Calculate posterior mean for a set of test data points for one region

INPUT:
 model - Considered probability model
 prediction - Vectors of realizations of the posterior distribution for the test data points
 param - Addiotnal parameter for the GPD model, otherwise not used
 number_policies - Number of trials if the considered model is of Binomial type
OUTPUT:
 Posterior mean for data points in vector 'expectation'
*/

void derive_expectation(const string &model, const vector<vector<double> > &prediction,
			const double &param, const vector<int> &number_trials,
			vector<double> &expectation){

  //Case 1: Binomial model with succes probability on logit scale
  if(model == "Binomial"){
    for(int i=0; i<prediction.size(); ++i){
      double temp = 0;
      for(int j=0; j<prediction[i].size(); ++j)
	temp = temp + exp(prediction[i][j])/(1+exp(prediction[i][j]));
      expectation.push_back(temp/(double)prediction[i].size() * (double)number_trials[i]);
    }
  }
  //Case 2:Poisson model with rate on log scale
  else if(model == "Poisson"){
    for(int i=0; i<prediction.size(); ++i){
      double temp = 0;
      for(int j=0; j<prediction[i].size(); ++j)
	temp = temp + exp(prediction[i][j]);
      expectation.push_back(temp/(double)prediction[i].size());
    }
  }
  //Case 3: Gaussian model
  else if(model == "Gaussian"){
    for(int i=0; i<prediction.size(); ++i){
      double temp = 0;
      for(int j=0; j<prediction[i].size(); ++j)
	temp = temp + prediction[i][j];
      expectation.push_back(temp/(double)prediction[i].size());
    }
  }
  //Case 4: truncated Binomial model with succes probability on logit scale
  else if(model=="tBinom"){
    for(int i=0; i<prediction.size(); ++i){
      double temp = 0;
      for(int j=0; j<prediction[i].size(); ++j){
	double temp2 = 1 - exp(-number_trials[i] * log(1 + exp(prediction[i][j])));
	temp = temp + exp(prediction[i][j])/(1+exp(prediction[i][j]))*1/temp2;
      }
      expectation.push_back(temp/(double)prediction[i].size() * (double)number_trials[i]);
    }
  }
  //Case 5: Generalized Pareto distribution
  else{
    for(int i=0; i<prediction.size(); ++i){
      double temp = 0;
      for(int j=0; j<prediction[i].size(); ++j)
	temp = temp + exp(prediction[i][j])/(1-param);
      expectation.push_back(temp/(double)prediction[i].size());
    }
  }
}
