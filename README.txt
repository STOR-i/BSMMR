This is the code for the Bayesian Spatial Monotonic Multiple Regression approach by Rohrbeck et al. and the .R used to create and analyse the results

Address for correspondence: 

Christian Rohrbeck
STOR-i CDT, Lancaster University, UK
c.rohrbeck@lancaster.ac.uk

------------------------------------------------------------------
In order to run the programme, the following steps are required
------------------------------------------------------------------

1) Compile the C++ files in "Code" using the provided Makefile in "Code" 

2) Copy the created executable "BSMMR" and the file "BayesOptim.R" into the folder where it is supposed to run

3) Create an empty folder called "Results"

4) Provide the required input files:
  (i) "input.txt" which specifies the parameters
  (ii) "data.txt" for the observations with the indices running from 1 to K
        (a) for "Gaussian", "Poisson", "GPD" the sequence per line is 
            Indice of the region - Observation y - Covariate observations x
        (b) for "Binomial" the sequence per line is
            Indice of the region - Observation y - Number of trials A - Covariate observations x
  (iii) "neighbours.txt" specifies the neighbourhood structure with each represententing one region. 
  (iv)  "weights.txt" specifies the weights d_{k,k'} in the prior density with each line referring to one region
  (v)   "weights_baseline.txt" contains the weigths for the CAR modelling of the baseline levels if necessary 
  (vi)  If available, one may also provide test data in a separate file called "test_data.txt" in the same format as "data.txt".
        The programm then returns the predictive mean squared error for the test data set.  
  
5) Run the executable "BSMMR" 

6) The estimated marked point processes and other parameters are stored in "Results" and posterior mean and convergence is provided in "Analysis"

------------------------------
Remarks on the current version
------------------------------
Compiling the provided could would result in a BSMMR implementation which ignores the baslines. In order to estimate them too,
the function "update_baseline" on line 125 in "main.cpp" needs to included. 
