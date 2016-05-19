#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <random>
#include "analysis.hpp"
#include "function_evaluation.hpp"

using namespace std;


/*Read MCMC output and call further functions to analyse the MCMC output and predict test data

INPUT:
 model                   - String describing the probability model
 dim                     - Dimension of the full covariate space
 K                       - Number of regions
 grid_size               - Grid points in each direction resulting in regular grid 
 DELTA_MIN               - Smallest possible level the function can take
 convergence_sample_size - Number of points sampled for each region to check convergence
 min_element             - Lower left corner of the full covariate spaces of the regions
 max_element             - Upper right corner of the full covariate spaces of the regions
 generator               - Random number generator             
OUTPUT:
 Analysed MCMC output in the folder Analysis 
SUBFUNCTIONS CALLED:
 check_convergence, function_on_grid, detect_thresholds, prediction
*/

void analysis_and_prediction(const string &model, const int&dim, const int &K,
			     const int &grid_size, const double &DELTA_MIN,
			     const int &convergence_samples_size, 
			     const vector<vector<vector<double> > > &sample_boundaries, 
			     mt19937 &generator){
 
  //Import test data
  vector<vector<vector<double> > > test_covariates(K, vector<vector<double> > ());
  vector<vector<int> > test_number_trials(K, vector<int> ());
  vector<vector<double> > true_response(K, vector<double> ());
  import_test_data(test_covariates, test_number_trials, true_response, dim, model);
  
  //Create folder "Analysis" and initialize file streams
  int temp_command = system("mkdir Analysis");
  ifstream baselinefile, processfile, parameterfile;
  ofstream processsizefile; 
  double temp; // Temporary variable to store values imported from .txt files

  //Read basline levels from "baseline.txt"
  vector<vector<double> > baseline( K, vector<double> () );
  vector<vector<double> > param( K, vector<double> () );
  parameterfile.open("Results/parameters.txt");
  baselinefile.open("Results/baseline.txt");
  while(parameterfile >> temp){
    for(int k=0; k<K; ++k){
      baselinefile >> temp;
      baseline[k].push_back(temp);
      if(model == "GPD" || model == "Gaussian"){
	parameterfile >> temp;
	param[k].push_back(temp);
      }
    }
  }

  //Read point process estimates from the different files
  int number_points;
  processsizefile.open("Analysis/Process_sizes.txt");
  
  //Perform analysis for each Region individually
  for(int k=0; k<K; ++k){
    cout<<"Region "<<k+1<<endl;

    vector<vector<vector<double> > > process_output;
    int pos = 0;
    string filename = "Results/Region" + to_string(k+1) + ".txt";

    //Read number of points in the process, store points and write the number to external file
    processfile.open(filename);
    while(processfile >> number_points){
      processsizefile << number_points << " ";
      process_output.push_back(vector<vector<double > > ());
      for(int i=0; i<number_points; ++i){ 
	vector<double> point;
	for(int j=0; j<dim+1; ++j){
	  processfile>>temp;
	  point.push_back(temp);
	}
	process_output[pos].push_back(point);
      }
      ++pos;
    }
    processsizefile<<"\n"; 
    processfile.close(); 

    //Check convergence of the chain
    check_convergence(convergence_samples_size, process_output, baseline[k], k, dim,
		      DELTA_MIN, sample_boundaries[k], generator);
   
    //Derive function on a regular grid of the full covariate space
    function_on_grid(process_output, baseline[k], dim, k, grid_size, DELTA_MIN,
		     sample_boundaries[k]);

    //Detect potential threshold effects
    detect_thresholds(process_output, baseline[k], sample_boundaries[k][0], dim, k,
		      DELTA_MIN);

    //Derive mean squared error of prediction
    prediction(process_output, baseline[k], param[k], test_covariates[k], test_number_trials[k],
	       true_response[k], dim, model, DELTA_MIN);
  }

  processsizefile.close();   
}

/*Check convergence by sampling a given number of points and derive assigned assigned levels
  for these points for each iteration step

INPUT:
 sample_size       - Number of points to be sampeld for checking convergence
 process           - Set of point processes describing the function for each iteration stored
 baseline          - Assigned baseline for each iteration step
 k                 - Indice of the region to be considered
 dim               - Dimension of the full covariate space
 DELTA_MIN         - Smallest possible level the point can take
 sample_boundaries - Lower left and upper right corner of the space on which function was estimated
 generator         - Random number generator
OUTPUT:
 MCMC chain for the sampled points are stored in "Analysis/Region(k)_convergence.txt"
SUBFUNCTIONS_CALLED:
 "function_evaluation.cpp":
  evaluate_test_point
*/

void check_convergence(const int &sample_size, const vector<vector<vector<double> > > &process_k,
		       const vector<double> &baseline_k, const int &k, const int &dim,
		       const double &DELTA_MIN, const vector<vector<double> > &sample_boundaries,
		       mt19937 &generator){

  //Open file stream to store results for the sample points
  ofstream convergencefile;
  convergencefile.open("Analysis/Region" + to_string(k+1) +"_convergence.txt"); 

  for(int i=0; i<sample_size; ++i){
    //Sample point for checking convergence uniformly across the sample space and store file
    vector<double> point;
    for(int j=0; j<dim; ++j){
      uniform_real_distribution<double> uniform(sample_boundaries[0][j], sample_boundaries[1][j]);
      point.push_back(uniform(generator));
      convergencefile<<point[j]<<" ";
    }

    //Store assigned level for the sampled point for each iteration in one line in the file
    for(int j=0; j<process_k.size(); ++j)
      convergencefile << evaluate_test_point(point, process_k[j], baseline_k[j], dim,
					     DELTA_MIN) <<" ";
    convergencefile<<"\n";
  }

  convergencefile.close();
}
 
/*Derive mean and quantiles of the estimated posterior function on a regular grid for Region k  

INPUT:
 process           - Set of point processes describing the function for each iteration stored
 baseline          - Assigned baseline level fore Region k for each iteration step
 dim               - Dimension of the full covariate space
 k                 - Indice of the region to be considered 
 grid_size         - Grid points in each direction resulting in regular grid
 DELTA_MIN         - Smallest possible level the point can take
 sample_boundaries - Lower left and upper right corner of the space on which function was estimated
OUTPUT:
 Derived posterior mena and quantiles are stored in "Region(k)_function.txt"
SUBFUNCTIONS CALLED:
 derive_summary_statistics, update_pos
*/

void function_on_grid(const vector<vector<vector<double> > > &process,
		      const vector<double> &baseline, const int &dim, const int &k,
		      const int &grid_size, const double &DELTA_MIN,
		      const vector<vector<double> > &sample_boundaries){

  //Open file to store results 
  ofstream gridfile;
  gridfile.open("Analysis/Region" + to_string(k+1)+"_function.txt");

  //Store header line in the external text file
  for(int j=0; j<dim; ++j)
    gridfile << "Covariate" + to_string(j+1) << " ";
  gridfile<< "Mean "<<"Median "<<"1%_quantile " <<"5%_quantile "<<"95%_quantile " <<
    "99%_quantile"<<"\n";

  //Derive step sizes in each direction
  vector<double> step_size (dim, 0);
  for(int i=0; i<dim; ++i)
    step_size[i] = (sample_boundaries[1][i] - sample_boundaries[0][i]) / grid_size;
  
  vector<int> pos(dim, 0);
  bool end_of_grid = false;
  while(end_of_grid != true){
    //Derive and store location of grid point
    vector<double> point;
    for(int i=0; i<dim; ++i){
      point.push_back( sample_boundaries[0][i] + pos[i]*step_size[i] );
      gridfile << point[i]<<" ";
    }

    //Store summary statistics for the grid point
    vector<double> result = derive_summary_statistics(point, process, baseline, dim, DELTA_MIN);
    for(int i=0; i<result.size(); ++i)
      gridfile << result[i]<<" ";
    gridfile<<"\n";

    //Derive next grid point, if not finished
    update_pos(pos, dim-1, grid_size, end_of_grid);
  }
  gridfile.close();
}


/*Derive next grid point to be considered by incrementing the indice in the jth component

INPUT:
 temp_pos    - Current grid cel considered
 j           - Indice of the covariate in which the next point lies
 grid_size   - Grid points in each direction resulting in regular grid with (grid_size)^dim points
 end_of_grid - Boolean variable to check whether all points have been considered
OUTPUT:
 Updated set of indices for the grid cell and updated bollean variable
*/

void update_pos(vector<int> &pos, const int &j, const int &grid_size, bool &end_of_grid){

  // Increment jth component by 1
  pos[j] = pos[j] + 1; 

  //Check whether end of grid in the jth direction is not reached 
  if(pos[j] <= grid_size) return;
  //If end of grid reached, update boolean variable to indicate that all grid cells are considered
  else if(pos[0] > grid_size) end_of_grid = true;
  //Otherwise consider the (j+1)th component
  else{
    pos[j] = 0;
    update_pos(pos, j-1, grid_size, end_of_grid);
  }
}

/*Derive mean and quantiles for a given grid point

INPUT:
 point      - Grid point in the full covariate space
 process_k  - Set of point processes describing the function for each iteration stored
 baseline_k - Assigned basleine level for Region k for each iteration step
 dim        - Dimension of the full covariate space
 DELTA_MIN  - Smallest possible level the point can take
OUTPUT:
 Sumamry statistics for the considered grid point
SUBFUNCTIONS CALLED:
 "function_evaluation.txt"
  evaluate_test_point
*/

vector<double> derive_summary_statistics(const vector<double> &point, 
					 const vector<vector<vector<double > > > &process_k, 
					 const vector<double> &baseline_k, const int &dim,
					 const double &DELTA_MIN){

 vector<double> result;

  //Derive assigned level at 'point'for each iteration stored
  vector<double> level;
  for(int i=0; i<process_k.size(); ++i){
    //If process is empty at one iteraion step, the level is equal to DELTA_MIN + baseline
    if(process_k[i].size() == 0) level.push_back(DELTA_MIN + baseline_k[i]);
    //Else derive the functional level
    else level.push_back(evaluate_test_point(point, process_k[i], baseline_k[i], dim, DELTA_MIN));
  }
 
  //Derive posterior mean
  double mean_level = 0;
  for(int i=0; i<level.size(); ++i)
    mean_level = mean_level + level[i];
  result.push_back( mean_level/level.size() );
  
  //Derive posterior median
  sort(level.begin(), level.end());
  if(level.size() %2 == 0)
    result.push_back( (level[level.size()/2 - 1] + level[level.size()/2])/2 );
  else
    result.push_back( level[(level.size() - 1)/2] );

  //Derive posterior 1%, 5%, 95% and 99% quantiles
  int pos_1 = level.size() * 0.01;
  result.push_back(level[pos_1]);
  int pos_5 = level.size() * 0.05;
  result.push_back(level[pos_5]);
  int pos_95 = level.size() * 0.95;
  result.push_back(level[pos_95]);
  int pos_99 = level.size() * 0.99;
  result.push_back(level[pos_99]);

  return result;
} 

/*Detection of potential threshold effects for Region k from the MCMC chains

INPUT:
 process_k   - Set of point processes describing the function for each iteration stored
 baseline    - Assigned baseline level for Region k for each iteration step
 min_element - Lower left corner of the sample space for region k
 dim         - Dimension of the full covariate space
 k           - Indice of the region to be considered 
 DELTA_MIN   - Smallest possible level the point can take
OUTPUT:
 Potential threshold effects and their probability of occurence are stored in external .txt file
SUBFUNCTIONS CALLED:
 evaluate_test_point
*/

void detect_thresholds(const vector<vector<vector<double> > > &process_k, 
		       const vector<double> &baseline_k, const vector<double> &min_element,
		       const int &dim, const int &k, const double &DELTA_MIN){

  vector<vector<double> > thresholds;
  for(int i=0; i<process_k.size(); ++i){
    for(int j=0; j<process_k[i].size(); ++j){
      
      //Shift point slightly such that previous location puts constraint on new point
      vector<double> point = process_k[i][j];
      int id = 0;
      while(0<1){
	if(point[id] >  min_element[id] + 1e-12){
	  point[id] = point[id] - 1e-12;
	  break;
	}
	else ++id;
      } 
      
      //Derive level of the shifted point
      double level_point = evaluate_test_point(point, process_k[i], baseline_k[i], dim, DELTA_MIN);
      
      //Check whether point is potential threshold and if so whether it is a new threshold or not
      if(process_k[i][j].back() + baseline_k[i] - level_point < 0.1) continue;
      else{//Otherwise, the point is a threshold and we check whether it's a new one or not
	bool temp = true;
	//Loop over thresholds found so far and check whether point is close to one of them
	for(int m=0; m<thresholds.size(); ++m){
	  //If level difference is too big, the points cannot describe the same threshold
	  if(abs(process_k[i][j].back() + baseline_k[i] - thresholds[m][dim]) > 0.01) continue;
	  else{//Otherwise, check whether they are close in the covariate space
	    double dist = 0.0;//Initialize distance
	    for(int n=0; n<dim; ++n)
	      dist  = dist + pow(process_k[i][j][n] - thresholds[m][n], 2);
	    //If distance too large, the point is a new threshold
	    if(sqrt(dist) > 0.02) continue;
	    //Otherwise, derive weighted average of the points in location and mark
	    else{
	      for(int n=0; n<dim; ++n)
		thresholds[m][n] = (thresholds[m].back()*thresholds[m][n]+process_k[i][j][n])/
		  (thresholds[m].back() + 1.0);	      
	      thresholds[m][dim] =
		(thresholds[m].back()*thresholds[m][dim] + process_k[i][j].back() + 
		 baseline_k[i])/(thresholds[m].back() + 1.0);
	      thresholds[m].back() = thresholds[m].back() + 1;
	      temp = false;
	      break;
	    }
	  }
	}
	//If point is new threshold effect, store it as new point
	if(temp == true){
	  vector<double> threshold_new = process_k[i][j];
	  threshold_new[dim] = threshold_new[dim] + baseline_k[i];
	  threshold_new.push_back(1.0);
	  thresholds.push_back(threshold_new);
	}
      }
    }
  }

  //Derive empirical occurrence probabilities 
  for(int i=0; i<thresholds.size(); ++i)
    thresholds[i].back() = thresholds[i].back()/process_k.size();

  //Export detected threshold effects and occurrence probability to file stream
  ofstream thresholdfile;
  thresholdfile.open("Analysis/Region" + to_string(k+1) + "_thresholds.txt");
  for(int j=0; j<dim; ++j){
    string temp_string = "Covariate" + to_string(j+1);
    thresholdfile << temp_string<<" ";
  }
  thresholdfile<< "Probability "<<"\n";
  for(int i=0; i<thresholds.size(); ++i){
    for(int j=0; j<thresholds[i].size(); ++j)
      thresholdfile << thresholds[i][j]<<" ";
    thresholdfile<<"\n";
  }
  thresholdfile.close();
}

/* Import test data set from external "test_data.txt" file

INPUT
 test_covariates    - Empty vector to be filled
 test_number_trials - Empty vector to be filled
 true_response      - Empty vector to be filled
 dim                - Dimension of the full covariate space
 model              - String describing the probability model
*/

void import_test_data(vector<vector<vector<double> > > &test_covariates,
		      vector<vector<int> > &test_number_trials,
		      vector<vector<double> > &true_response,
		      const int &dim, const string &model){

  int k, number_trials;
  double response;
  vector<double> covariate( dim, 0.0 );

  ifstream datafile;
  datafile.open("Input/test_data.txt");
  // Each line gives the response, the number of trials and covariate values for a test data point
  while(datafile >> k){
    datafile >> response;
    true_response[k-1].push_back(response);
    if(model == "Binomial" || model == "tBinom" || model =="Bernoulli"){
      datafile >> number_trials;
      test_number_trials[k-1].push_back(number_trials);
    }
    for(int i=0; i<dim; ++i)
      datafile >> covariate[i];
    test_covariates[k-1].push_back(covariate);
  }
  datafile.close();
}

/* Derive mean predictive error for a set of test points

INPUT:
 process            - Set of processes representing the function at each step of the Markov chain
 baseline           - Assigned baseline levels at each of the setps 
 param              - Additional stream of model parameters (if reqired for probability model)
 test_covariates    - Covariate observations for the test points
 test_number_trials - Number of trials for the tes points (for Binomial, tBinom or Bernoulli)
 true_response      - Observations related to the covariate observations in the test data set
 dim                - Dimension of the full covariate space
 model              - String describing the probability model
 DELTA_MIN          - Lowest possible level the point can take
OUTPUT:
 Mean squared error of predicted mean and true observation of test data points printed in terminal
SUBFUNCTIONS CALLED:
 evaluate_test_point, derive_predictive_mean
*/

void prediction(const vector<vector<vector<double> > > &process, const vector<double> &baseline,
		const vector<double> &param, const vector<vector<double> > &test_covariates,
		const vector<int> &test_number_trials, const vector<double> &true_response,
		const int &dim, const string &model, const double &DELTA_MIN){

  double MSE;
  vector<double> predictive_mean;
  vector<vector<double> > levels( true_response.size(), vector<double> (process.size(),0.0) );

  // Derive assigned level for each of the test data points at each step of the MCMC chain
  for(int i=0; i<true_response.size(); ++i){
    for(int j=0; j<process.size(); ++j)
      levels[i][j] =
	evaluate_test_point(test_covariates[i], process[j], baseline[j], dim, DELTA_MIN);
  }

  // Derive predictive mean for all test data points
  derive_predictive_mean(model, levels, param, test_number_trials, predictive_mean);

  // Calculate mean squared error of true response and predicted mean and print it in Terminal
  for(int i=0; i<true_response.size(); ++i)
    MSE = MSE + pow( predictive_mean[i] - true_response[i] , 2.0 );
  cout <<"Predictive squared error: "<< MSE <<" for "<< true_response.size()<<" points."<< endl;  
}

/*Derive predictive mean for multiple points based on their assigned levels from the MCMC
  algorithm

INPUT:
 model           - String describing the probability model
 levels          - Assigned chains of functional levels for each of the test data points
 param           - Additional stream of model parameters in the probability model (if required)  
 number_trials   - Number of trials (for Binomial, tBinomial or Bernoulli random variables
 predictive_mean - Empty vector to be filled
OUTPUT
 Predictive mean for each point returned in the 'predictive_mean' vector
*/

void derive_predictive_mean(const string &model, const vector<vector<double> > &levels,
			    const vector<double> &param, const vector<int> &number_trials,
			    vector<double> &predictive_mean){

  //Case 1: Binomial model with succes probability on logit scale
  if(model == "Binomial"){
    for(int i=0; i<levels.size(); ++i){
      double prob = 0;
      for(int j=0; j<levels[i].size(); ++j)
	prob = prob + 1/( 1 + exp(-levels[i][j]) );
      predictive_mean.push_back(prob / (double)levels[i].size() * (double)number_trials[i]);
    }
  }
  //Case 2:Poisson model with rate on log scale
  else if(model == "Poisson"){
    for(int i=0; i<levels.size(); ++i){
      double mean = 0;
      for(int j=0; j<levels[i].size(); ++j)
	mean = mean + exp(levels[i][j]);
      predictive_mean.push_back(mean / (double)levels[i].size());
    }
  }
  //Case 3: Gaussian model
  else if(model == "Gaussian"){
    for(int i=0; i<levels.size(); ++i){
      double mean = 0;
      for(int j=0; j<levels[i].size(); ++j)
	mean = mean + levels[i][j];
      predictive_mean.push_back(mean / (double)levels[i].size());
    }
  }
  //Case 4: truncated Binomial model with succes probability on logit scale
  else if(model=="tBinom"){
    for(int i=0; i<levels.size(); ++i){
      double mean = 0;
      for(int j=0; j<levels[i].size(); ++j){
	double prob_excess = 1 - exp(-number_trials[i] * log(1 + exp(levels[i][j])) );
	mean = mean + exp( levels[i][j])/(1 + exp(levels[i][j]) ) * 1/prob_excess;
      }
      predictive_mean.push_back(mean / (double)levels[i].size() * (double)number_trials[i]);
    }
  }
  //Case 5: Generalized Pareto distribution
  else if(model=="GPD"){
    for(int i=0; i<levels.size(); ++i){
      double mean = 0;
      for(int j=0; j<levels[i].size(); ++j)
	mean = mean + exp(levels[i][j]) / (1 - param[j]);
      predictive_mean.push_back(mean / (double)levels[i].size());
    }
  }
  //Case 6: Bernoulli random variable 
  else{
    for(int i=0; i<levels.size(); ++i){
      double mean = 0.0;
      for(int j=0; j<levels[i].size(); ++j){
	mean = mean + 1 - exp( -(double)number_trials[i] * log( 1 + exp(levels[i][j])) );
      }
      predictive_mean.push_back(mean / (double)levels[i].size());
    }
  }
}

