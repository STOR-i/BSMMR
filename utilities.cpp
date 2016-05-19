#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>

using namespace std;

/*Read the parameters for RJMCMC algorithm, model and analysis*/

void read_configurationfile(int &number_iterations, int &burn_in, int &size_thinning,
			    double &p_birth, double &p_death, int &max_jumps, string &model,
			    double &alpha_param, double &beta_param, int &K, int &dim,
			    double &DELTA_MIN, double &DELTA_MAX, string &set, string &cross_val,
			    double &omega, double &omega_0, int &number_repetitions,
			    double &p, double &q, double &eta, double &tau_0, double &alpha_tau,
			    double &beta_tau, int &grid_size, int &convergence_samples_size){
 
  ifstream configurationfile; 
  configurationfile.open("Input/input.txt"); 
  configurationfile.ignore(256, ':');
  configurationfile>>number_iterations;
  configurationfile.ignore(256, ':');
  configurationfile>>burn_in;
  configurationfile.ignore(256, ':');
  configurationfile>>size_thinning;
  configurationfile.ignore(256, ':');
  configurationfile>>p_birth;
  configurationfile>>p_death;
  configurationfile.ignore(256, ':');
  configurationfile>>max_jumps;
  configurationfile.ignore(256, ':');
  configurationfile>>model;
  if(model=="Gaussian"|| model=="GPD"){
    configurationfile>>alpha_param;
    configurationfile>>beta_param;
  }
  configurationfile.ignore(256, ':');
  configurationfile>>K;
  configurationfile.ignore(256, ':');
  configurationfile>>dim;
  configurationfile.ignore(256, ':');
  configurationfile>>DELTA_MIN;
  configurationfile>>DELTA_MAX;
  configurationfile.ignore(256, ':');
  configurationfile>>set;
  configurationfile.ignore(256, ':');
  configurationfile>>cross_val;
  if(cross_val=="false") configurationfile>>omega;
  else{
    configurationfile>>omega_0;
    configurationfile>>number_repetitions;
  }
  configurationfile.ignore(256, ':');
  configurationfile>>p;
  configurationfile>>q;
  configurationfile.ignore(256, ':');
  configurationfile>>eta;
  configurationfile.ignore(256, ':');
  configurationfile>>tau_0;
  configurationfile>>alpha_tau;
  configurationfile>>beta_tau;
  configurationfile.ignore(256, ':');
  configurationfile>>grid_size;
  configurationfile.ignore(256, ':');
  configurationfile>>convergence_samples_size;
  configurationfile.close();
}

/*Import training data to fit the model and initialise variables associated to the data

INPUT:
 covariates          - Empty matrix to be filled with covariate observations for K regions
 reponse             - Empty matrix to be filled with observations of the response variable
 number_trials       - Empty matrix to be filled with number of trials if model is "Binomial"
 neighbours          - Empty matrix to be filled with the neighbourhood structure
 weights             - Empty matrix matrix to be filled with the weights d_k,k' in prior density
 weights_baseline    - Empty matrix to be filled with the weights for smoothing baseline levels 
 assigned_level_curr - Empty matrix which containing current levels for the covariate observations
 dim                 - Dimension of the full covariate space
 DELTA_MIN           - Smallest possible level
 model               - String defining the probability model
OUTPUT:
 Updated matrices
*/

void import_data(vector<vector<vector<double> > > &covariates, vector<vector<double> > &response,
		 vector<vector<int> > &number_trials, vector<vector<int> > &neighbours,
		 vector<vector<double> > &assigned_level_curr, vector<vector<double> > &weights,
		 vector<vector<double> > &weights_baseline, const int &dim,
		 const double &DELTA_MIN, const string &model, const int &K){		 
 
  covariates = vector<vector<vector<double> > > (K, vector<vector<double> > () );
  response = vector<vector<double> > (K, vector<double> () );
  number_trials = vector<vector<int> > (K, vector<int> () );

  int region, temp_number_trials;
  double temp_response;
  vector<double> covariate(dim, 0.0);

  //Read data observations from "datafile.txt"
  ifstream datafile;
  datafile.open("Input/data.txt");
  //Read indice of the region the observation comes from and attach values to respective entry
  while(datafile >> region){
    datafile>>temp_response;
    response[region-1].push_back(temp_response);
    if(model == "Binomial" || model=="tBinom" || model=="Bernoulli"){
      datafile>>temp_number_trials;
      number_trials[region-1].push_back(temp_number_trials);
    }
    for(int i=0; i<dim; ++i)  datafile>>covariate[i];
    covariates[region-1].push_back(covariate);
  }
  datafile.close();

  //As functions are constant to DELTA_MIN at the start, assign each data point level DELTA_MIN
  for(int k=0; k<K; ++k)
    assigned_level_curr.push_back( vector<double> (response[k].size(), DELTA_MIN) );

  //Read neighbourhood structure and weigths for the spatial smoothing
  ifstream neighbourdata, weightsdata,  weightsbaselinedata;
  neighbourdata.open("Input/neighbours.txt");
  weightsdata.open("Input/weights.txt");
  weightsbaselinedata.open("Input/weights_baseline.txt");

  int size;
  //Read number of neighbours from "neighbours.txt" and then indices of the neighbours and weights
  while(neighbourdata >> size){
    vector<int> neighbours_k(size, 0);
    vector<double> weight1(size, 0);
    vector<double> weight2(size, 0);
    for(int i=0; i<size; ++i){
      neighbourdata>>neighbours_k[i];
      --neighbours_k[i];
      weightsdata>>weight1[i];
      weightsbaselinedata>>weight2[i];
    }
    neighbours.push_back(neighbours_k);
    weights.push_back(weight1);
    weights_baseline.push_back(weight2);
  }
  neighbourdata.close();
  weightsdata.close();
  weightsbaselinedata.close();
}

/*Initialize the parameters to run the RJMCMC algorithm*/

void initialise_parameters(vector<vector<vector<vector<double> > > > &processes,
			   vector<vector<int> > &ind, vector<double> &baseline,
			   vector<double> &param, vector<vector<int> > &low_sets,
			   vector<vector<int> > &upp_sets, vector<vector<vector<int> > > &pos_cov,
			   vector<vector<double> > &sum_levels,
			   vector<vector<int> > &number_jumps_subprocesses,
			   vector<int> &number_jumps,vector<double> &mean,
			   vector<double> &sum_weights_baseline,  
			   const vector<vector<double> > &weights_baseline, const int &dim,
			   const int &K){

  //Initialize indices of the covariate subsets - ind
  for(int i=0; i<dim; ++i){
    vector<vector<int> > temp_ind = ind;
    for(int j=0; j<ind.size(); ++j){
      vector<int> temp =ind[j];
      temp.push_back(i);
      temp_ind.push_back(temp);
    }
    ind = temp_ind;
    ind.insert(ind.begin() + i, vector<int>({i}));
  }
  
  //Initialize variables specifying point processes
  baseline = vector<double>(K, 0.0);
  param = vector<double>(K, 0.1);
  number_jumps = vector<int> (K, 0);
  mean = vector<double> (K, 0.0);
  for(int k=0; k<K; ++k){
    processes.push_back(vector<vector<vector<double> > > (ind.size(), vector<vector<double> > ()));
    sum_levels.push_back(vector<double>(ind.size(), 0) );
    number_jumps_subprocesses.push_back(vector<int> (ind.size(), 0));
  }
  
  //Initialize low_sets, upp_sets and pos_cov for comparing points from different subsets
  low_sets = vector<vector<int> > (ind.size(), vector<int> () );
  upp_sets = vector<vector<int> > (ind.size(), vector<int> () );
  pos_cov = vector<vector<vector<int> > > (ind.size(), vector<vector<int> > ()) ;
  for(int j=0; j<ind.size(); ++j) pos_cov[j] = vector<vector<int> > (ind.size(), vector<int> ());
  
  for(int i=0; i<dim; ++i){
    low_sets[i].push_back( {i} );
    upp_sets[i].push_back( {i} );
    for(int j=0; j<ind.size(); ++j){
      if(i==j) pos_cov[i][j].push_back( {0} );
      else pos_cov[j][i].push_back( {-1} );
    }
  }

  for(int i=dim; i<ind.size(); ++i){
    for(int j=0; j<ind.size(); ++j){
      bool temp=true;
      for(int k=0; k<ind[j].size(); ++k){
	if(find(ind[i].begin(), ind[i].end(), ind[j][k]) != ind[i].end()){
	  pos_cov[j][i].push_back(find(ind[i].begin(), ind[i].end(), ind[j][k]) - ind[i].begin());
	  continue;
	}
	else{
	  temp=false;
	  break;
	}
      }
      if(temp == true){
	low_sets[i].push_back( {j} );
	upp_sets[j].push_back( {i} );
      }
      else{
	pos_cov[j][i] = vector<int>( {-1} );
      }
    }
  }

  // Initialize sum of the weights for each region
  for(int k=0; k<K; ++k){
    double temp =0;
    for(int i=0; i<weights_baseline[k].size(); ++i)
      temp = temp + weights_baseline[k][i];
    sum_weights_baseline.push_back(temp);
  }
}

/*Derive lower left and upper right corner of the spaces on which the proposed new points are
  sampled and derive the associated voume of the sample space in dependency on the considered
  domain of integration  

INPUT:
 covariates - Covariate observations for the K regions
 dim        - Dimension of the full cvariate space
 K          - Number of regions
 set        - String describing the domain on which the difference of the functions is evaluated
 ind        - Index set describing the considered covariate subsets
 neighbours - Neighbourhood structure of the regions
OUTPUT:
 A                   - Lower left and upper right corner of A_k,k' for each pair of neighbours k_k'
 sample_boundary     - Lower left and upper right conrer of sample space for each region k
 volume_sample_space - Volume of sample space with respect to covariate subset on which proposed
*/

void derive_volume(vector<vector<vector<vector<double> > > > &A,
		   vector<vector<vector<double> > > &sample_boundaries,
		   vector<vector<double> > &volume_sample_space,
		   const vector<vector<vector<double> > > &covariates,
		   const int &dim, const int &K, const string &set,
		   const vector<vector<int> > &neighbours, const vector<vector<int> > &ind){

  A = vector<vector<vector<vector<double> > > > (K, vector<vector<vector<double> > > ());
  sample_boundaries = vector<vector<vector<double> > > (K, vector<vector<double> > ());

  //Derive minimum and maximum covariate observation in each covariate for each region 
  vector<vector<vector<double> > > bound_covariates(K, vector<vector<double> > () );
  for(int k=0; k<K; ++k){
    bound_covariates[k] = vector<vector<double> > (2, covariates[k][0] );
    for(int i=1; i<covariates[k].size(); ++i){
      for(int j=0; j<dim; ++j){
	bound_covariates[k][0][j] = min(covariates[k][i][j], bound_covariates[k][0][j]);
	bound_covariates[k][1][j] = max(covariates[k][i][j], bound_covariates[k][1][j]);
      }
    }
  }

  //Derive sample_boundaries and A in dependency on whether and how covariate space are varying  

  //Case 1: The covariate spaces are set to be equal
  if(set=="Equal"){
    //Derive bounds in each covariate separately
    vector<double> temp_min = bound_covariates[0][0];
    vector<double> temp_max = bound_covariates[0][1];    
    for(int k=1; k<K; ++k){
      for(int j=0; j<dim; ++j){
	temp_min[j] = min(temp_min[j], bound_covariates[k][0][j]);
	temp_max[j] = max(temp_max[j], bound_covariates[k][1][j]);	
      }
    }
    sample_boundaries =
      vector<vector<vector<double> > > (K, vector<vector<double> > ( {temp_min, temp_max} ) );
    for(int k=0; k<K; ++k)  A[k].push_back(sample_boundaries[0]);
  }

  //Case 2: The integral is considered on the union of the covariate spaces X_k and X_k'
  else if(set=="Union"){
    for(int k=0; k<K; ++k){
      sample_boundaries[k]=
	vector<vector<double> > ( { bound_covariates[k][0], bound_covariates[k][1] } );
      //Derive A_k,k' for each neighbour k' of Region k
      for(int j=0; j<neighbours[k].size(); ++j){
	vector<double> temp_min = bound_covariates[k][0];
	vector<double> temp_max = bound_covariates[k][1];
	//Derive bounds in each covariate separately
	for(int m=0; m<dim; ++m){
	  temp_min[m] = min(temp_min[m], bound_covariates[neighbours[k][j]][0][m]);
	  temp_max[m] = max(temp_max[m], bound_covariates[neighbours[k][j]][1][m]);
	  sample_boundaries[k][0][m] =
	    min(sample_boundaries[k][0][m], bound_covariates[neighbours[k][j]][0][m]);
	  sample_boundaries[k][1][m] =
	    max(sample_boundaries[k][1][m], bound_covariates[neighbours[k][j]][1][m]);
	}
	A[k].push_back(vector<vector<double> > ( {temp_min, temp_max} ) );
      }
    }
  }

  //Case 3: The integral is evaluated only over the overlap of the covariate spaces X_k and X_k'
  else{
    for(int k=0; k<K; ++k){
      sample_boundaries[k] =
	vector<vector<double> > ( { bound_covariates[k][0], bound_covariates[k][1] } );
      //Derive A_k,k' for each neighbours k' of Region k
      for(int j=0; j<neighbours[k].size(); ++j){
	vector<double> temp_min = bound_covariates[k][0];
	vector<double> temp_max = bound_covariates[k][1];
	//Derive bounds in each covariate separately
	for(int m=0; m<dim; ++m){
	  temp_min[m] = max(temp_min[m], bound_covariates[neighbours[k][j]][0][m]);
	  temp_max[m] = min(temp_max[m], bound_covariates[neighbours[k][j]][1][m]);
	}
	A[k].push_back(vector<vector<double> > ( {temp_min, temp_max} ) );
      }
    }
  }

  //Derive the volume of the sample space with respect to the subspace the proposal is made on
  for(int k=0; k<K; ++k){
    volume_sample_space.push_back(vector<double> ());
    for(int j=0; j<ind.size(); ++j){
      double temp = 1.0; 
      for(int m=0; m<ind[j].size(); ++m){
	temp = temp * (sample_boundaries[k][1][ind[j][m]] - sample_boundaries[k][0][ind[j][m]]);
      }
      volume_sample_space[k].push_back(temp);
    }
  }
}

/*Store the current processes, random effects and smoothing parameter to different files

INPUT:
 rep           - Number of the current repetition
 processes     - Marked point processes for all K regions
 K             - Number of regions
 dim -         - Dimension of thw full covariate space
 ind           - Index set describing the considered covariate subsets
 random_effect - Vector of random effects assigned to the K regions
 tau           - Spatial smoothing parameter of the random effects
 number_jumps  - Vector containing the number of points describing the function for each region
 min_elment    - Lopwer left corner of the sample spaces of the point processes
*/

void write_to_file(const vector<vector<vector<vector<double> > > > &processes, const int &K,
		   const int &dim, const vector<vector<int> > &ind,
		   const vector<double> &baseline, const double &tau, const vector<double> &param,
		   const string &model, const vector<int> &number_jumps,
		   const vector<vector<vector<double> > > &sample_boundaries){

  //Store samples in separate file for each region 
  for(int k=0; k<K; ++k){
    string filename = "Results/Region" + to_string(k+1) + ".txt";
    ofstream process_file;
    process_file.open(filename, ios::app);
    process_file<<number_jumps[k] << "\n";//Store number of points describing monotonic function
   
    //Translate each point to lie in the full covairate space and attach it to the file
    for(int i=0; i<processes[k].size(); ++i){
      for(int j=0; j<processes[k][i].size(); ++j){
	vector<double> point_full = sample_boundaries[k][0];
	for(int m=0; m<ind[i].size(); ++m)
	  point_full[ind[i][m]] = processes[k][i][j][m];
	for(int m=0; m<dim; ++m)	  
	  process_file << point_full[m] << " ";	
	process_file << processes[k][i][j].back() << "\n";
      }
    }
    process_file.close();
  }
 
  //Store baseline for all regions in one file
  ofstream baseline_file;
  baseline_file.open("Results/baseline.txt", ios::app);
  for(int k=0; k<K; ++k)
    baseline_file << baseline[k]<<" ";
  baseline_file << "\n";
  baseline_file.close();

  //Store remaining updated parameters in additional .txt file
  ofstream parameter_file;
  parameter_file.open("Results/parameters.txt", ios::app);
  parameter_file << tau << " ";//Store the current spatial smoothing parameter of the baseline
  if(model=="Gaussian" || model=="GPD"){
    for(int k=0; k<K; ++k)
      parameter_file << param[k] << " ";
  }
  parameter_file << "\n";
  parameter_file.close();
}
