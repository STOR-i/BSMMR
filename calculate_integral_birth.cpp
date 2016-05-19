#include <iostream>
#include <vector> 
#include <algorithm>
#include <cmath>
#include "calculate_integral_birth.hpp"
#include "function_evaluation.hpp"

using namespace std;

/*Calculate prior ratio for proposing a 'Birth' for Region k
  
INPUT:
 set               - Indicating whether covariate spaces are varying or not
 k                 - Indice of the region which is updated
 proposal          - Proposed location and functional level
 ip                - Indice of the subset on which the proposal is made
 processes         - Point processes describing the monotonic functions for all regions
 ind               - Indices of the covariate subsets on which the processes are defined 
 neighbours        - Vector containing the indices of the neighbours of Region k 
 weights           - Constants d_k,k' specifying similarity of Region k and its neighbours
 dim               - Dimension of the full covariate space
 DELTA_MIN         - Smallest possible level
 sample_boundaries - Lower and upper corner on which points are proposed
 A_k               - Left and upper corner of A_k,k' for all neighbours k' of Region k
 low_sets          - Vector containing the subprocesses on which the proposal may put constraint
 pos_cov           - Gives the indices to compare if points are defined on different subsets
 level_curr        - Current level at the proposed location
 p,q               - Values for the distance measure
OUTPUT: 
 Integrated difference of the proposed process with the neighbours -
 integrated difference of the current process with the neighbours   
SUBFUNCTIONS CALLED:
  simplify_process, derive_current_process, integral_difference 
 "function_evaluation.hpp":
  function_evaluation, evaluate_at_proposal, check_monotonicity, monotonicity_full,
*/

double integral_birth(const string &set, const int &k, const vector<double> &proposal,
		      const int &ip, const vector<vector<vector<vector<double> > > > &processes,
		      const vector<vector<int> > &ind, const vector<int> &neighbours,
		      const vector<double> &weights, const int &dim, const double &DELTA_MIN,
		      const vector<vector<double> > &sample_boundaries,
		      const vector<vector<vector<double> > > &A_k, const vector<int> &low_sets,
		      const vector<vector<vector<int> > > &pos_cov, const double &level_curr,
		      const double &p, const double &q){

  double integral = 0.0;

  //Translate proposal to lie in the full covariate space with current level as mark
  vector<double> proposal_full = sample_boundaries[0];
  for(int m=0; m<ind[ip].size(); ++m) proposal_full[ind[ip][m]] = proposal[m];
  proposal_full.push_back( level_curr );

  /*If covariate spaces are not varying, the proposal has not to be translated further,
    leading to higher efficiency in the integral calculation*/
  if(set=="Equal"){
    
    //Derive minimum set of points giving area affected by the proposal and the current function
    vector<vector<double> > border, temp_border, process_k;
    vector<double> bound = A_k[0][1];
    simplify_process (proposal_full, proposal, processes[k], ind, dim, level_curr, pos_cov[ip],
		      process_k, border, temp_border, bound);
  
    //Derive functions over the area affected by the proposal for the neighbours of Region k
    vector<vector<vector<double> > > processes_neigh;
    for(int i=0; i<neighbours.size(); ++i){
      double temp_level = level_at_proposal(processes[neighbours[i]], proposal, ip, DELTA_MIN,
					    low_sets, pos_cov);
      processes_neigh.push_back( vector<vector<double> > () );
      derive_process_neigh(proposal_full, proposal, processes[neighbours[i]], ind, dim, temp_level,
			   pos_cov[ip], border, bound, temp_border, processes_neigh[i]);
    }
    //Calculate ratio in the prior densities based on the derived point processes
    integral_difference(process_k, processes_neigh, border, temp_border, dim, proposal.back(),
			ind.back(), DELTA_MIN, p, q, weights, integral);
  }

  /*If covariate spaces are varying, proposal has to be translated individually for each pair*/
  else{
    for(int i=0; i<neighbours.size(); ++i){
      int temp_ip =ip;
      double temp_level_curr = level_curr;
      vector<double> temp_proposal_full = proposal_full;
      vector<double> temp_proposal = proposal;

      //Check whether proposal affects area A_k,k' over which the integral is evaluated
      //If upper corner of A_k,k' puts no constraint on proposed location, integral won't change
      if( check_monotonicity(proposal, A_k[i][1], ind[ip], 0) == false ) continue;      
      //Else check whether proposal has to be translated
      if( monotonicity_full(A_k[0][1], proposal_full, dim) == false ){
	//If its ha to be translated, check whether proposal effects A_k,k'
	for(int m=0; m<dim; ++m) temp_proposal_full[m] = max( A_k[i][0][m], proposal_full[m] );
	temp_level_curr = function_evaluation(processes[k], ind, dim, temp_proposal_full,
					      level_curr);
	if(temp_level_curr > proposal.back()) continue;
	//Update variables for calculating integral difference, if proposal affects A_k,k' 
	else{
	  temp_proposal_full.back() = temp_level_curr;
	  temp_ip = ind.size() - 1;
	  temp_proposal = temp_proposal_full;
	  temp_proposal.back() = proposal.back();
	}
      }
      
      //Derive minimum set of points giving area affected by the proposal and the current function
      vector<double> bound = A_k[i][1];
      vector<vector<double> > border, temp_border, process_k;
      simplify_process (temp_proposal_full, temp_proposal, processes[k], ind, dim, temp_level_curr,
			pos_cov[temp_ip], process_k, border, temp_border, bound);

      //Derive functions over the area affected by the proposal for the neighbour i of Region k
      vector<vector<vector<double> > > processes_neigh;
      processes_neigh.push_back(vector<vector<double> > ());      
      double temp_level = function_evaluation(processes[neighbours[i]], ind, dim,
					      temp_proposal_full, DELTA_MIN);
      derive_process_neigh(temp_proposal_full, temp_proposal, processes[neighbours[i]], ind, dim,
			   temp_level,pos_cov[temp_ip], border, bound, temp_border,
			   processes_neigh[0]);
     
      //Calculate ratio in the prior densities based on the derived point processes
      integral_difference(process_k, processes_neigh, border, temp_border, dim,
			  temp_proposal.back(), ind.back(), DELTA_MIN, p, q, weights, integral);
    }
  }
  return(integral);
}

/*Derive corner points of the area affected by the proposal and the current function over it

INPUT: 
 proposal_full - Translated proposal which lies in the full covariate space
 proposal      - Proposed location and mark
 process       - Set of point processes describing the monotonic function
 ind           - Indices describing the covariate sets considered
 dim           - Dimension of the full covariate space
 level_curr    - Current level of the monotonic function at the location of the proposal
 pos_cov       - Indices used to compare points which are defined on different covariate subsets
 process_k     - Empty point process
 border        - Empty point process
 bound         - Upper corner of the covariate space
OUTPUT:
 process_k - Point process describing the current function over the region affected by the proposal
 border    - Point process desribing the area affected by the proposal
 bound     - Highest covariate values over all border points
SUBFUNCTIONS CALLED:
 check_monotonicity, monotonicity_sets ("function_evaluation.hpp")
 minimum_border, minimum_process_check_border, sort_by_level
*/

void simplify_process(vector<double> &proposal_full, const vector<double> &proposal,
		      const vector<vector<vector<double> > > &process,
		      const vector<vector<int> > &ind, const int &dim, const double &level_curr, 
		      const vector<vector<int> > &pos_cov, vector<vector<double> > &process_k,
		      vector<vector<double> > &border, vector<vector<double> > &temp_border,
		      vector<double> &bound){ 

  //Set value at proposal to current value and attach point to process describing current function
  process_k.push_back( proposal_full );

  //STEP 1: Consider point processes defined on one covariate only
  for(int i=0; i<dim; ++i){
    vector<double> temp_point = proposal_full;
    for(int j=0; j<process[i].size(); ++j){
      //While mark smaller than currentl level at the proposal, ignore the point
      if(process[i][j].back() <= level_curr) continue;
      //If mark higher than proposed level, update upper bound and stop considering process
      else if( process[i][j].back() > proposal.back() ){
	bound[i] = process[i][j][0];
	break;
      }
      //Otherwise, store points as candidates for describing current process over that area
      else{
	temp_point[i] = process[i][j][0];
	temp_point.back() = process[i][j].back();
	process_k.push_back( temp_point );
      }
    }
    //Store border point with updated upper bound
    temp_point[i] = bound[i];
    border.push_back( temp_point );
  }

  //STEP 2: Consider point processes defined on more than one covariate
  for(int i=dim; i<process.size(); ++i){
    vector<double> temp_point = proposal_full;
    //If points in the process may put monotonic constraint on proposal, make additional check
    if(pos_cov[i][0] != -1){
      for(int j=0; j<process[i].size(); ++j){
	//If mark smaller than current level, condsider next point 
	if( process[i][j].back() <= level_curr ) continue;
	//Check whether point lies within the boundaries 
	else if( process[i][j][0] > bound[ind[i][0]] ) break;
	else if( check_monotonicity(process[i][j], bound, ind[i], 1) == false ) continue;
	/*Check whether point puts constraint on 'proposal' and store it as border point
	  in a temporary process*/
	else if( monotonicity_sets(proposal, process[i][j], 0, pos_cov[i]) == true){
	  for(int m=0; m<ind[i].size(); ++m)
	    temp_point[ind[i][m]] = process[i][j][m];
	  temp_border.push_back( temp_point );
	}
	/*Else translate point to lie on the edge of the area affected by the proposal and
	  store it in either 'process_k' or 'border' dependent on its mark */
	else{
	  for(int m=0; m<ind[i].size(); ++m)
	    temp_point[ind[i][m]] = max( proposal_full[ind[i][m]], process[i][j][m] );
	  if( process[i][j].back() > proposal.back() )
	    border.push_back( temp_point );
	  else{
	    temp_point.back() = process[i][j].back();
	    process_k.push_back( temp_point );
	  }   
	}
      }
    }

    //Else potential points can only lie on the edge of the area affected by the proposal
    else{
      for(int j=0; j<process[i].size(); ++j){
	//If mark smaller than current level, condsider next point 
	if( process[i][j].back() <= level_curr ) continue;
	//Check whether point lies within the boundaries 
	else if( process[i][j][0] > bound[ind[i][0]] ) break;
	else if( check_monotonicity(process[i][j], bound, ind[i], 1) == false ) continue;
      	/*Else translate point to lie on the edge of the area affected by the proposal and
	  store it in either 'process_k' or 'border' dependent on its mark */
	else{
	  for(int m=0; m<ind[i].size(); ++m)
	    temp_point[ind[i][m]] = max( proposal_full[ind[i][m]], process[i][j][m] );
	  if( process[i][j].back() > proposal.back() )
	    border.push_back( temp_point );
	  else{
	    temp_point.back() = process[i][j].back();
	    process_k.push_back( temp_point );
	  }   
	}
      }
    }
  }
 
  //STEP 3: Derive border process and point process decribing the current function over that area
  //Derive minimum point process describing the area affected by the proposal
  minimum_border(border, temp_border, bound, dim);
  //Delete points in the current process which do not describe the current function
  minimum_process_check_border(process_k, border, dim);
  //Sort the points describing the current function by level from smallest to highest
  sort_by_level(process_k);
}	       		       					    

/*Derive set of points describing area affected by the proposal from set of potential candidates

INPUT
 border      - Point process with original points putting no constraint on translated  proposal
 temp_border - Point process with originial points putting constraint on translated proposal
 bound       - Maximum covariate values border point can take
 dim         - Dimension of the full covariate space
OUTPUT
 Minimized processes border and temp_border such that no border point puts monotonic constraint
 on another border point and updated 'bound'
SUBFUNCTIONS CALLED
 monotonicity_full("function_evaluation.hpp") 
*/

void minimum_border(vector<vector<double> > &border, vector<vector<double> > &temp_border,
		    vector<double> &bound, const int &dim){

  //STEP 1: Find minimum values for 'bound' and remove points which are outside of these
  for(int i=0; i<dim; ++i){ 
    for(int j=dim; j<border.size(); ++j){
      //If component i of point j larger than highest possible vovariate value, delete point j 
      if(border[i][i] < border[j][i]){
	border.erase( border.begin() + j );
	--j;
      }
      //If smaller value for bound i is found, replace point i by point j
      else if( monotonicity_full(border[j], border[i], dim) == true ){
	border[i][i] = border[j][i];
	border.erase( border.begin() + j );
	--j;
      }
      else continue;
    }
  }
  for(int i=0; i<dim; ++i) bound[i] = border[i][i];

  //STEP 2: Check monotonic constraints for remaining points
  for(int i=dim; i<border.size(); ++i){
    for(int j=0; j<i; ++j){
      //If point j puts constraint on point i, delete point j 
      if( monotonicity_full(border[i], border[j], dim) == true ){
	border.erase( border.begin() + j );
	--j;
	--i;
      }
      //If point i puts constraint on point j, delete point i 
      else if( monotonicity_full(border[j], border[i], dim) == true ){
	border.erase( border.begin() + i );
	--i;  
	break;
      }
      else continue;
    }
  }

  //STEP 3: Consider the points which put monotonic constraint on the proposal
  for(int i=0; i<border.size(); ++i){ 
    for(int j=0; j<temp_border.size(); ++j){
      //If a point from temp_border puts constraint on a point in border, delete the former
      if( monotonicity_full(border[i], temp_border[j], dim) == true ){
	temp_border.erase( temp_border.begin() + j );
	--j;
      }
      else continue;
    }
  }

  //STEP 4:Check whether point in temp_border puts constraint on another point in temp_border
  for(int i=1; i<temp_border.size() ; ++i){
    //Check whether two points in temp_border put monotonic constraint on each other
    for(int j=0; j<i; ++j){
      //If point j puts constraint on point i, delete point j 
      if( monotonicity_full(temp_border[i], temp_border[j], dim) == true ){
	temp_border.erase( temp_border.begin() + j );
	--i;
      }
      //If point i puts constraint on point j, delete point i 
      else if( monotonicity_full(temp_border[j], temp_border[i], dim) == true ){
	temp_border.erase( temp_border.begin() + i );
	--i;  
	break;
      }
      else continue;
    }
  }
}

/*Derive set of points describing current function over the area affected by the proposal

INPUT
 process - Set of poential points describing the current function
 border  - Set of points descrbing the area affected by the proposal
 dim     - Dimension of the full covariate space
OUTPUT
 Set of process describing the current function over the are described by 'border'
SUBFUNCTIONS CALLED
 minimum_process, check_border
*/

void minimum_process_check_border(vector<vector<double> > &process, 
				  const vector<vector<double> > &border, const int &dim){

  //Delete points which do not describe current function
  minimum_process(process, dim, 2);

  //Check whether points lie within the area described by border and the proposal
  for(int i=1; i<process.size(); ++i){
    //If a candiate puts constraint on a border point, it lies outside the area and is deleted
    if( check_border(process[i], border, dim, 0) == false){
      process.erase( process.begin() + i );
      --i;
    }
    else continue;
  }
}

/*Delete points which are not describing the current function over the area affected the proposal

INPUT
 process - Set of points which are potentially describing the current function
 dim     - Dimension of the full covariate space
OUTPUT
 Cleared process describing the current function
SUBFUNCTIONS CALLED
 monotonicity_full ("function_evaluation.hpp")
*/

void minimum_process(vector<vector<double> > &process, const int &dim, const int &pos){

  for(int i=pos; i<process.size(); ++i){
    for(int j=0; j<i; ++j){
      //If point j puts constraint on point i put the mark of point i is higher, delete point j 
      if( monotonicity_full(process[i], process[j], dim) == true &&
	  process[i].back() >= process[j].back() ){
	process.erase( process.begin() + j);
	--i;
	--j;
      }
     //If point i puts constraint on point j put the mark of point i is higher, delete point i 
      else if( monotonicity_full(process[j], process[i], dim) == true &&
	      process[j].back() >= process[i].back()){
	process.erase( process.begin() + i );
	--i;
	break;
      }
      else continue;
    }
  }
} 

/*Check whether a candidate descrbining the current function over the area affected by the
  proposal lies outside the area by checking monotonicity for all but the first pos border points

INPUT
 point  - Candidate for describing the current function
 border - Points describing the area ffecetd by the proposal
 dim    - Dimension of the full covariate space
 pos    - Indice from which border point comparison is started
 */

bool check_border( const vector<double> &point, const vector<vector<double> > &border,
		   const int &dim, const int &pos){

  for(int i=pos; i<border.size(); ++i){
    //If the point puts monotonic constraint on border point i, return 'false'
    if( monotonicity_full(border[i], point, dim) == true) return false;
    else continue;//Consider next border point
  }
  return true;//If no monotonic constraint for any border point, return 'true'
}


/*Sort points by their mark

INPUT
 process - Point process with last component giving the mark
OUTPUT
 Sort process from smallest to highest mark
*/

void sort_by_level(vector<vector<double> > &process){

  for(int i=1; i<process.size()-1; ++i){
    //Find indice of point with ith smallest level
    int pos = i;
    for(int j=i+1; j<process.size(); ++j){
      if( process[pos].back() > process[j].back() ) pos = j;
      else continue;
    }
    //Swap point such that point with ith lowest level is a position i after iteration i
    process[i].swap(process[pos]);
  }
}

/*Derive function over an area given the border points and the proposal for neighbour 

INPUT:
 process   - Set of processes which describes the monotonic function
 ind - Indices describing the covariate sets considered
 proposal_transformed - proposal transformed to lie in the full covariate space
 border - process containing the set of points desribing the area affected by the proposal
 max_element - right upper corner of the full covarite space

 pos_cov    - Matrix used to conpare points from two different set. Gives the indices to compare.
 proposal  - Original proposed location and mark
 ind_prop  - Indice of the subset the proposal is made on
 process_current - empty process
OUTPUT:
 process_current - Set of points describing function over the region affected by the proposal
SUBFUNCTIONS CALLED:
  level_at_proposal, check_monotonicity, check_monotonicty_subsets
*/

void derive_process_neigh(vector<double> &proposal_full, const vector<double> &proposal,
			  const vector<vector<vector<double> > > &process,
			  const vector<vector<int> > &ind, const int &dim,const double &temp_level,
			  const vector<vector<int> > &pos_cov,
			  const vector<vector<double> > &border, const vector<double> &bound,
			  const vector<vector<double> > &temp_border, 
			  vector<vector<double> > &process_k2){

  //Set mark at proposal equal to current level and attach point to current process
  proposal_full.back() = temp_level;
  process_k2.push_back( proposal_full );

  //Initalize temporary process
  vector<vector<double> > temp_process;

  //Derive set of potential points describing the function over the area affceted by the proposal
  for(int i=0; i<process.size(); ++i){
    vector<double> temp_point = proposal_full;
    //If points in the process may put monotonic constraint on proposal, make additional check
    if(pos_cov[i][0]!=-1){
      for(int j=0; j<process[i].size(); ++j){
	//If mark smaller than current level, consider next point
	if(process[i][j].back() <= temp_level) continue;
	//Check whether point lies within the boundaries
	else if(process[i][j][0] > bound[ind[i][0]]) break;
	else if( check_monotonicity(process[i][j], bound, ind[i], 1) == false ) continue;
	//Check whether point puts constraint on 'proposal'
	else if( monotonicity_sets(proposal, process[i][j], 0, pos_cov[i]) == true ){	  
	  //If so, translate point and check border condition for border points in 'temp_border'
	  for(int m=0; m<ind[i].size(); ++m)
	    temp_point[ind[i][m]] = process[i][j][m];
	  if( check_border(temp_point, temp_border, dim, 0) == true ){
	    if( check_border(temp_point, border, dim, dim) == true ){
	      temp_point.back() = process[i][j].back();
	      process_k2.push_back( temp_point );	    
	    }
	    else continue;
	  }
	  else continue;
	} 
	//Else translate point and check border condition for points in 'border'
	else{
	  for(int m=0; m<ind[i].size(); ++m) 
	    temp_point[ind[i][m]] = max( proposal_full[ind[i][m]], process[i][j][m] );
	  if( check_border(temp_point, border, dim, dim) == true ){
	    temp_point.back() = process[i][j].back();
	    temp_process.push_back( temp_point );
	  }
	  else continue;
	}
      }
    }

    //Else potential points can only lie on the edge of the area affected by the proposal
    else{
      for(int j=0; j<process[i].size(); ++j){
	//If mark smaller than current level, condsider next point 
	if(process[i][j].back() <= temp_level) continue;
	//Check whether point lies within the boundaries 
	else if(process[i][j][0] > bound[ind[i][0]]) break;
	else if( check_monotonicity(process[i][j], bound, ind[i], 1) == false ) continue;
	/*Else translate point to lie on the edge of the area affected by the proposal and
	  check border condition for points in 'border' */
	else{
	  for(int m=0; m<ind[i].size(); ++m) 
	    temp_point[ind[i][m]] = max( proposal_full[ind[i][m]], process[i][j][m] );
	  if( check_border(temp_point, border, dim, dim) == true ){
	    temp_point.back() = process[i][j].back();
	    temp_process.push_back( temp_point );
	  }
	}
      }
    }
  }

  //STEP 2: Remove points violating monotonic constraint from the process
  minimum_process(temp_process, dim, 0);
  //Merge the two processes 
  process_k2.insert( process_k2.end(), temp_process.begin(), temp_process.end() );
  //Sort points with increasing mark
  sort_by_level(process_k2);
}

/*Calculate the integrated difference of the process the proposal and its neighbours by
  considering each point in the original process seprately.

INPUT:
 process_k   - Process describing the current function on the area the proposal is made on
 processes_neigh - Processes of the neighbouring functions over the area affected by proposal
 border             - Set of points describing the area affacted by the proposal
 proposal           - Level of the proposal
 dim                - Dimension of the full covariate space
 ind                - Indices of the covariates
 DELTA_MIN          - Smallest possible level
 p,q                - Values for the distance measure
 weights            - Constants describing the prior belief of the two processes being similar
 integral           - Initial vale of the integral
OUTPUT:
 Derived integral difference for the acceptance probability
SUBFUNCTIONS CALLED:
 check_monotonicty, function_evaluation_volume, calculate_volume
*/

void integral_difference(vector<vector<double> > &process_k, 
			 vector<vector<vector<double> > > &processes_neigh,
			 vector<vector<double> > &border, vector<vector<double> > &temp_border_old,
			 const int &dim, const double &proposal, vector<int> ind, 
			 const double &DELTA_MIN, const double &p, const double &q,
			 const vector<double> &weights, double &integral){ 

  //Merge the two border processes
  border.insert( border.begin(), temp_border_old.begin(), temp_border_old.end() );

  //Consider each point in the process describing the current function separately
  for(int i=process_k.size()-1; i>=0; --i){
    //Derive area over the border area which has level equal to point process_k[i]
    vector<vector<double> > temp_border;
    derive_temp_border(process_k[i], border, dim, ind, process_k[0], temp_border);

    //Consider the neighbours separately as prior is defined pair-wise
    for(int k=0; k<processes_neigh.size(); ++k){
      
      //Derive function for the neighbouring regions over the area temp_border
      vector<vector<double> > temp_process;
      //Copy border for next neighbour
      vector<vector<double> > temp_border2 = temp_border;
      
      //Derive level at process_k[i] and store point as the first one in temp_process
      double temp_level = evaluate_at_point(process_k[i], processes_neigh[k], dim);
      vector<double> temp_point = process_k[i];
      temp_point.back() = temp_level;
      temp_process.push_back(temp_point);
      //Loop over the remaining points and derive process over the area given bu temp_border
      for(int j=processes_neigh[k].size() - 1; j>0; --j){
	//If level of point smaller than temp_level, ignore point
	if(processes_neigh[k][j].back() <= temp_level) break;
	//If point puts monotonic constraint on process_k[i], move it to temporary process
	else if( monotonicity_full(process_k[i], processes_neigh[k][j],dim) == true ){
	  temp_process.insert(temp_process.begin() + 1, processes_neigh[k][j]);
	  processes_neigh[k].erase( processes_neigh[k].begin() + j );
	}
	/*Else translate the point such that it lies on the edge of the area affected by
	  process_k[i] and check border condition*/
	else{
	  for(int m=0; m<dim; ++m)
	    temp_point[m] = max( process_k[i][m], processes_neigh[k][j][m] );
	  if( check_border(temp_point, temp_border, dim, 0) == true ){ 
	    temp_point.back() = processes_neigh[k][j].back();	    
	    temp_process.insert( temp_process.begin() + 1, temp_point );
	  }
	  else continue;
	}
      }
      //Remove redundant points from the temporary process
      minimum_process(temp_process, dim, 1);

      //Calculate the integral difference for each level in temp_process
      for(int j=temp_process.size() - 1; j>=0; --j){ 
	//Derive area with level of process_k[i] and temp_process[j]
	vector<vector<double> > temp_border3;
	for(int m=0; m<temp_border2.size(); ++m){
	  //If border point puts constraint on temp_process[j], move it to the temporary border
	  if( monotonicity_full(temp_process[j], temp_border2[m], dim) == true ){
	    temp_border3.push_back( temp_border2[m] );
	    temp_border2.erase( temp_border2.begin() + m );
	    --m;
	  }
	  else{
	  //Else translate it to lie on the edge of the area of interest and 
	    vector<double> temp_point = temp_process[j];
	    for(int n=0; n<dim; ++n) 
	      temp_point[n] = max(temp_border2[m][n], temp_process[j][n]);
	    temp_border3.push_back(temp_point);
	  }
	}
	//Remove redundant points and attach temp_process[j] to both border processes
	minimum_region(temp_border3, dim);
	temp_border3.insert( temp_border3.begin(), temp_process[j] );
	temp_border2.insert( temp_border2.begin(), temp_process[j] );

	//Derive integral difference for area with the levels of process_k[i] and temp_process[j]
	integral = integral + weights[k] * calculate_volume(temp_border3, dim) *
	  (abs(pow(1+temp_process[j][dim]-DELTA_MIN, p)-pow(1+proposal-DELTA_MIN, p))* 
	   pow(abs(temp_process[j][dim] - proposal), q) - 
	   abs(pow(1+temp_process[j][dim]-DELTA_MIN, p)-pow(1+process_k[i][dim]-DELTA_MIN, p))* 
	   pow(abs(temp_process[j][dim] - process_k[i][dim]), q));
      }
    }
    //After evaluating whole area with level equal to process_k[i], attach point to border process
    border.push_back(process_k[i]);
  }
}

/*Derive set of border points describing the area with level equal to level at a point

INPUT
 point       - Point representing current function over a part of area affected by the proposal
 border      - Point process describing the area affected by the proposal
 dim         - Dimension of the full covariate space
 ind         - Indices for the full covariate space
 proposal    - Proposed new level
 temp_border - Empty point process
OUTPUT
 Minimum point process describing the area with level equal to the one of point in 'temp_border'
SUBFUNCTIONS CALLED
 monotonicity_full ("function_evaluation.cpp")
*/

void derive_temp_border(const vector<double> &point,vector<vector<double> > &border,
			const int &dim, vector<int> &ind, vector<double> &proposal,
			vector<vector<double> > &temp_border){

  //Initalize two additional  temporary border processes
  vector<vector<double> > temp_border1, temp_border2;

  //Detect covariate for which 'point' has same covariate value as 'proposal' and delete i from ind
  int pos = 0;
  while( point[pos] > proposal[pos] ) ++pos;
  ind.erase( ind.begin() + pos );

  //Derive border points and assign them to one of the border processes
  for(int i=0; i<border.size(); ++i){
    //If border point puts constraint on point,move it to the temporary border process
    if( monotonicity_full(point, border[i], dim) == true ){ 
      temp_border.push_back( border[i] );
      border.erase( border.begin() + i );
      --i;
    }
    /*Assign point to border  dependent on whether the original border point has the
      same covariate value than 'point' at component i*/
    else{
      vector<double> temp_point = border[i];
      for(int j=0; j<ind.size(); ++j)
	temp_point[ind[j]] = max( point[ind[j]], border[i][ind[j]] );
      if(border[i][pos] > point[pos]) temp_border1.push_back( temp_point );
      else temp_border2.push_back( temp_point );
    }
  }

  //Derive minimum set of points describing the area with level equal to the one of point
  minimum_region(temp_border1, dim);
  minimum_region(temp_border2, dim);
  //Reinsert indice deleted and merge all three border processes
  ind.insert(ind.begin() + pos, pos);
  temp_border.insert( temp_border.end(), temp_border1.begin(), temp_border1.end() );
  temp_border.insert( temp_border.end(), temp_border2.begin(), temp_border2.end() );
}
 
/*Derive maximum set of points such that no point puts monotonic constraint on another point

INTPUT
 border - Set of potential points for describing the area of interest
 dim    - Dimension of the full covariate space
OUTPUT
 Minimized process
SUBFUNCTIONS CALLED
 monotonicity_full ("functiona_evaluation.cpp")
*/
   
void minimum_region(vector<vector<double> > &border, const int &dim){ 

  //Consider each pair of points
  for(int i=1; i<border.size(); ++i){
    for(int j=0; j<i; ++j){
      //If point j puts constraint on point i, remove point j
      if( monotonicity_full(border[i], border[j], dim) == true ){
	border.erase( border.begin() + j );
	--j;
	--i;
      }
      //If point i puts constraint on point j, remove point i
      else if( monotonicity_full(border[j], border[i], dim) == true ){
	border.erase( border.begin() + i );
	--i;
	break;
      }
      else continue;
    }
  }
}    

/* Calculation of the volume of a multi-dimensional polytope with given border points by
1. Check whether the polytope is one-dimensional and derive smallest non-zero difference 
2. If polytope is two-dimensional, the volumn is derived directly
3. Consider the multi-dimensional case by:
 (a) Delete redundant border points, i.e.\ points which put a monotonic constraint on anther point
 (b) Sort the points by their first covariate
 (c) Slice the polytope based on this ordering and call function recursively for lower dimension

INPUT:
 process_border - Points of the process describing the polytope with first point being the proposal
 dim            - Dimension of the polytope
OUPUT: 
 volume - Volume of the polytope
SUBFUNCTIONS CALLED:
 monotonicty_full ("function_evaluation.cpp")
*/

double calculate_volume(vector<vector<double> > &border, int dim){
  
  double volume = 0.0;

  //If one-dimensional space, return volume directly
  if(dim == 1) volume = border[1][0] - border[0][0];

  //If space is two-dimensional
  else if(dim == 2){
    //Order points with respect to the first component
    for(int i=1; i<border.size()-1; ++i){
      int pos = i;
      for(int j=i+1; j<border.size(); ++j){
	if(border[pos][0] > border[j][0]) pos = j;
	else continue;
      }
      border[i].swap( border[pos] );
    }
    //Using the derived ordering, the volume is derived directly
    for(int i=2;  i<border.size(); ++i)
      volume = volume + (border[i-1][1] - border[0][1]) * (border[i][0] - border[i-1][0]);
  }

  //Otherwise, the polytope is defined in higher-dimensions
  else{
    //Order points with respect to the first component 
    for(int i=1; i<border.size()-1; ++i){
      int pos = i;
      for(int j=i+1; j<border.size(); ++j){ 
	if(border[pos][0] > border[j][0]) pos = j;
	else continue;
      }
      border[i].swap( border[pos] );
    }
    
    //Start by putting all point with smallest first component in one process
    int pos = 0;
    vector<vector<double> > temp_border;
    while( border[pos][0] == border[0][0] ){
      temp_border.push_back( vector<double> (border[pos].begin() + 1, border[pos].end()) );
      ++pos;
    }
    
    //Then add the next point with then smallest first component to the process 
    for(int i=pos; i<border.size(); ++i){
      
      /*As the second to dim th component stay the same, we can simplify the calculation to a 
	problem with only dim-1 dimensions*/
      volume = volume + (border[i][0] - border[i-1][0]) * calculate_volume(temp_border, dim-1);

      //After attaching new border point, remove redundant points from the border process
      temp_border.push_back( vector<double> (border[i].begin() + 1, border[i].end()) );      
      for(int j=1; j<temp_border.size()-1; ++j){
	if( monotonicity_full(temp_border.back(), temp_border[j], dim-1) == true ){
	  temp_border.erase( temp_border.begin() + j );
	  --j;
	}
	else continue;
      }
    } 
  }
  return volume;
}

