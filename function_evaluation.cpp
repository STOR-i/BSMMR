#include <iostream>
#include <vector>
#include <algorithm>
#include "function_evaluation.hpp"

using namespace std;

/*Derve functional level of the monotonic function at point arg

INPUT
 process - Set of point processes ordered by first component giving function
 ind     - Index set for the point processes
 dim     - Dimension of the full covariate space
 arg     - Point at which function to be evaluated
 level   - Lowest possible level
OUTPUT
 result  - Derived functional level for point arg
SUBFUNCTIONS CALLED:
 check_monotonicity
*/

double function_evaluation(const vector<vector<vector<double> > > &process,
			   const vector<vector<int> > &ind, const int &dim,
			   const vector<double> &arg, const double &level){

  double result = level;//Initalize return value
  
  //STEP 1: Consider one-dimensional subspaces
  for(int i=0; i<dim; ++i){ 
    int pos = 0;
    //Loop over process until covariate value is larger than the ith component of arg[i]
    for(int j=0; j<process[i].size(); ++j){
      if(process[i][j][0] < arg[i]) ++pos;
      else break;
    }
    //Check whether a point with smaller component was found and update 'result'
    if(process[i].size() > 0 && pos > 0)
      if(process[i][pos-1].back() > result)
	result = process[i][pos-1].back();
  }

  //STEP 2: Consider point processes defined on more than one covariate
  for(int i=dim; i<process.size(); ++i){
    for(int j=0; j<process[i].size(); ++j){
      //If level smaller than current level, consider next point       
      if(process[i][j].back() <= result) continue;
      //If first component larger than respective one of arg, consider next process
      else if(process[i][j][0] > arg[ind[i][0]]) break;
      //If 'arg' puts monotonic constraint on point, update level 
      else if( check_monotonicity(process[i][j], arg, ind[i], 1) == true )
	result = process[i][j].back();
      else continue;
    }
  }
  return result;
}

/*Derive functional levels for a set of points 

INPUT
 process   - Point processes describing function
 ind       - Covariate sets on which processes are defined on
 dim       - Dimension of the full covariate space
 args      - Points in the full covariate space for which levels are to be derived
 DELTA_MIN - Lowest possible functional level
OUTPUT
 Derived functional levels for points in 'args'
SUBFUNCTIONS CALLED
 check_monotonoicity
*/

vector<double> function_evaluation_set(const vector<vector<vector<double> > > &process,
				       const vector<vector<int> > &ind, const int &dim,
				       const vector<vector<double> > &args,
				       const double &DELTA_MIN){

  vector<double> result (args.size(), DELTA_MIN);//Initialize return variable

  //Consider each point separately
  for(int i=0; i<args.size(); ++i){
 
    //STEP 1: Consider point processes defined on one covariate only
    for(int j=0; j<dim; ++j){
      int pos = 0;
      //Loop over process until covariate value is larger than the ith component of arg[i]
      for(int k=0; k<process[j].size(); ++k){
	if(process[j][k][0] < args[i][j]) ++pos;
	else break;
      }
      //Check whether a point with smaller component was found and update 'result[i]'
      if(process[j].size() > 0 && pos > 0)
	if(process[j][pos-1].back() > result[i])
	  result[i] = process[j][pos-1].back();
    }
    
    //STEP 2: Consider point processes defined on more than one covariate
    for(int j=dim; j<process.size(); ++j){
      for(int k=0; k<process[j].size(); ++k){
	//If level smaller than current level, consider next point       
	if(process[j][k].back() <= result[i]) continue;
	//If first component larger than respective one of arg, consider next process
	else if(process[j][k][0] > args[i][ind[j][0]]) break;
	//If 'arg' puts monotonic constraint on point, update level 
	else if( check_monotonicity(process[j][k], args[i], ind[j], 1) == true )
	  result[i] = process[j][k].back();
	else continue;
      }
    }
  }
  return result;
}

/*Check whether monontonic criterion is fulfilled for all but the first pos components

INPUT:
 point - Point in subprocess with index set 'ind' 
 arg   - Point in the full covariate space
 ind   - Indices of the covariate subspace
 pos   - Component from which comparison is started
OUTPUT:
 'true' if arg puts monotonic constraint on point and 'false' otherwise
*/

bool check_monotonicity(const vector<double> &point, const vector<double> &arg,
			const vector<int> &ind, const int &pos){

  //Check for all but the first pos components monotonicity component-wise
  for(int i=pos; i<ind.size(); ++i){
    if(point[i] <= arg[ind[i]]) continue;//Monotonicity true in component, next  
    else return false;//Otherwise, monotonicity is violated and 'false' is returned
  }
  return true;//If all components fulfill monontonic criterion, return 'true'
}
 
/* Derive the minimum and maximum for sampling a proposed level at sampled location 

INPUT
 process    - Set of point processes describing the current process 
 proposal   - Location of the proposal
 ip         - Indice of the covariate subset on which proposal is made
 DELTA_MIN  - Minimum possible level
 DELTA_MAX  - Maximum possible level
 low        - Index sets on which 'proposal' may put monotonic constraint 
 upp        - Index sets which may put monotonic constraint on 'proposal' 
 pos_cov    - Matrix giving indices to compare point from different subspaces
OUTPUT
 result - Lower bound and upper bound for the level at point 'proposal' 
SUBFUNCTIONS CALLED: 
 monotonicty_sets
*/

vector<double> bounds_level(const vector<vector<vector<double> > > &process,
			    const vector<double> &proposal, const int &ip,
			    const double &DELTA_MIN, const double &DELTA_MAX,
			    const vector<int> &low, const vector<int> &upp,
			    const vector<vector<vector<int> > > &pos_cov){

  double lower_bound = DELTA_MIN; // Initalize lower bound for level
  double upper_bound = DELTA_MAX; // Initalize upper bound for level

  //Loop over the processes on which points in process 'ip' may put constraint
  for(int i=0; i<low.size(); ++i){ 
    for(int j=0; j<process[low[i]].size(); ++j){
      //If functional level smaller than current, consider next point
      if(process[low[i]][j].back() <= lower_bound) continue;
      //If first component greater than respective component of 'proposal' consider next process
      else if(process[low[i]][j][0] > proposal[pos_cov[low[i]][ip][0]]) break;
      //If proposal puts constraint on point j in process low[i], update lower bound
      else if(monotonicity_sets(process[low[i]][j], proposal, 1, pos_cov[low[i]][ip]) == true)
	lower_bound = process[low[i]][j].back();
      else continue;
    }
  }

  //Loop over the processes containing points which may put constraint on proposal
  for(int i=0; i<upp.size(); ++i){
    for(int j=0; j<process[upp[i]].size(); ++j){
      //If level of point higher than current upper bound, consider next point
      if(process[upp[i]][j].back() >= upper_bound) continue;
      //If point j in process upp[i] puts constraint on proposal, update upper bound 
      else if(monotonicity_sets(proposal, process[upp[i]][j], 0, pos_cov[ip][upp[i]]) == true)
	upper_bound = process[upp[i]][j].back();
      else continue;
    }
  }

  return vector<double> ( {lower_bound, upper_bound} );
}

/*Check monontonic constraint for two points from potentially different subprocesses

INPUT
  point1  - Point in a subprocess
  point2  - Point which is at least defined on the same components as 'point1' 
  pos     - Component from which to start comparison of components
  pos_cov - Matrix giving indices which have to be compared
OUTPUT
 'true' if 'point2' puts a monotonic constraint on 'point1' and 'false' otherwise
*/

bool monotonicity_sets(const vector<double> &point1, const vector<double> &point2,
		       const int &pos, const vector<int>& pos_cov){

  for(int i=pos; i<pos_cov.size(); ++i){//Compare all but the first pos components
    if(point1[i]<=point2[pos_cov[i]]) continue;//Check monotonic constraint 
    else return false;// If monotonicity is violated for one covariate, return 'false'
  }
  return true;//Return 'true', if condition is satisfied for all components
}

/* Derive the current value of the function at a point defined on a covariate subset 

INPUT
 process   - Point processes describing current function
 proposal  - Location of proposal
 ip        - Indices on which proposal is made
 DELTA_MIN - Minimum possible level
 low       - Index sets on which 'proposal' may put monotonic constraint 
 pos_cov   - Matrix giving indices to compare point from different subspaces
OUTPUT
 result    - Current level at the proposal
SUBFUNCTIONS CALLED 
 monotonicty_sets
*/

double level_at_proposal(const vector<vector<vector<double> > > &process,
			 const vector<double> &proposal, const int &ip, const double &DELTA_MIN,
			 const vector<int> &low, const vector<vector<vector<int> > > &pos_cov){

  double level = DELTA_MIN;//Initalize return variable

  //Loop over the processes on which points in process 'ip' may put constraint
  for(int i=0; i<low.size(); ++i){ 
    for(int j=0; j<process[low[i]].size(); ++j){//Loop over the points in process
      //If functional level smaller than current, consider next point
      if(process[low[i]][j].back() <= level) continue;
      //If first component greater than respective component of 'proposal' consider next process
      else if(process[low[i]][j][0] > proposal[pos_cov[low[i]][ip][0]]) break;
      //Check monotonicity of the point and the proposal on the subset defined by 'ind[low[i]]'
      else if(monotonicity_sets(process[low[i]][j], proposal, 1, pos_cov[low[i]][ip]) == true)
	level = process[low[i]][j].back();
      else continue;
    }
  }

  return level;
}

/* Derive sample boundaries for new location for proposed switch. Boundaries are
   derived such that the monotonicity structure is maintained.

INPUT
 ip          - Indice of the subset on which a shift is proposed
 point       - Point in process 'ip' to be shifted
 process     - Point processes describing the current function
 ind         - Index set associated with the covariate subsets
 min_element - Left lower corner in the sample space 
 max_element - Right upper corner in the sample space
 low         - Index sets on which 'proposal' may put monotonic constraint 
 upp         - Index sets which may put monotonic constraint on 'proposal' 
 pos_cov    - Matrix giving indices to compare point from different subspaces
OUTPUT
  boundaries   - Boundaries for sampling new location for point 'idp' 
*/

vector<vector<double> > bounds_location(const int &ip, const vector<double> &point,
					const vector<vector<vector<double> > > &process,
					const vector<vector<int> > &ind, 
					const vector<vector<double> > &sample_boundaries,
					const vector<int> &low, const vector<int> &upp,
					const vector<vector<vector<int> > > &pos_cov){

  //Initialize bounds and set them to smallest and highest possible points
  vector<vector<double> > bounds;
  for(int i=0; i<ind[ip].size(); ++i)
    bounds.push_back( { sample_boundaries[0][ind[ip][i]], sample_boundaries[1][ind[ip][i]] } );
  
  for(int i=0; i<low.size(); ++i){//Loop over the point processes
    for(int k=0; k<ind[low[i]].size(); ++k){//Loop over the components
      for(int j=0; j<process[low[i]].size(); ++j){//Loop over the points
	//Check whether point provides new lower bound in component
	if(process[low[i]][j][k] > bounds[ pos_cov[low[i]][ip][k] ][0] &&
	   process[low[i]][j][k] < point[ pos_cov[low[i]][ip][k] ])  
	  bounds[pos_cov[low[i]][ip][k]][0]= process[low[i]][j][k];
	//Check whether point provides new upper bound in component
	else if(process[low[i]][j][k]< bounds[ pos_cov[low[i]][ip][k] ][1] &&
		process[low[i]][j][k] > point[ pos_cov[low[i]][ip][k] ])  
	  bounds[pos_cov[low[i]][ip][k]][1]= process[low[i]][j][k];
	else continue;//Consider next point
      }
    }
  }
  
  //Consider point processes which may put monotonic constraint on proposal
  for(int i=1; i<upp.size(); ++i){//Loop over the point processes
    for(int k=0; k<ind[ip].size(); ++k){//Loop over the components
      for(int j=0; j<process[upp[i]].size(); ++j){//Loop over the points
	//Check whether point provides new lower bound in component
	if(process[upp[i]][j][pos_cov[ip][upp[i]][k]] > bounds[k][0] &&
	   process[upp[i]][j][pos_cov[ip][upp[i]][k]] < point[k]) 
	  bounds[k][0] = process[upp[i]][j][pos_cov[ip][upp[i]][k]];	
	//Check whether point provides new upper bound in component
	else if(process[upp[i]][j][ pos_cov[ip][upp[i]][k] ] < bounds[k][1] &&
		process[upp[i]][j][ pos_cov[ip][upp[i]][k] ] > point[k]) 
	  bounds[k][1] = process[upp[i]][j][pos_cov[ip][upp[i]][k]];
	else continue;//Consider next point
      }
    }
  }
  return bounds;//Return derived boundaries on the location
}


/*Derive level at a given point for a process described on the full covariate space only

INPUT:
 point   - Point in the full covariate space
 process - Process on the full covariate space (levels are ordered from small to high)
 dim     - Dimension of the full covariate space
OUTPUT:
 Level of the function described by the process at the point 
SUBFUNCTIONS CALLED:
 monotonicity_full
*/

double evaluate_at_point(const vector<double> &point, const vector<vector<double> > &process,
			 const int &dim){

  for(int j=process.size()-1; j>0; --j){//Loop over the points, starting with the highest level
    if( monotonicity_full(process[j], point, dim) == true )//Check monotonicity
      //If condition fulfilled, return level of the point, as all following have smaller level    
      return process[j].back();
  } 
  //If none of the 2nd to last highest point meet criterion, return smallest possible level
  return process[0].back();
}

/*Check whether point2 puts constraint on Point 1 with both defined on full covariate space

INPUT:
 point1 - Point 1 in full covariate space
 point2 - Point 2 in full covariate space
 dim    - Dimension of the full covariate space
OUTPUT
 'true' if Point 1 is smaller or equal than Point 2 in each covariate and 'false' otherwise
*/

bool monotonicity_full(const vector<double> &point1, const vector<double> &point2,
		       const int &dim){
  for(int i=0; i<dim; ++i){//Loop over all covariates
    if(point1[i] <= point2[i]) continue;//If monotonicity fulfilled, consider next component
    else return false;//If condition is violated, return 'false'
  }
  return true;//If monotonicity true for all covariates, return 'true'
}

/*Derive level for a point with marked point process defined on full covariate space

INPUT:
 point         - Point in the full covariate space to be evaluated
 process       - Marked point process describing monotonic function
 random_effect - Assigned random effect in the model
 dim           - Dimension of the full covariate space
 DELTA_MIN     - Smallest possible level the point can take
OUTPUT: 
 Derived sum of 'baseline' and level of the function at 'point' 
SUBFUNCTIONS CALLED:
 monotonicity_full
*/

double evaluate_test_point(const vector<double> &point, const vector<vector<double> > &process,
			   const double &baseline, const int &dim, const double &DELTA_MIN){

  double level = DELTA_MIN;

  for(int i=0; i<process.size(); ++i){
    //If level of point smaller than current level, consider next point
    if(process[i].back() <= level) continue;
    //Otherwise, check monotonic constraint and update assigned level if necessary
    else if(monotonicity_full(process[i], point, dim) == true) level = process[i].back();
    else continue;
  }
  return level + baseline;//Return derived level plus baseline level
}

