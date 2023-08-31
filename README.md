# Bayesian Sample Size Determination for Multilevel Models with Longitudinal Data

In this project, multilevel data under two hypotheses is simulated to which a multilevel model is fitted. Bayes factors for each hypotheses are calculated and their dependence on sample size (and other quantities of interest) is investigated. The hypothesis $H_0$ claims that there is no effect of the intervention (i.e., no slope difference in the two groups) while $H_1$ postulates that people in the intervention condition get better faster. In terms of model parameters, this means that $H_0$ states that the interaction parameter between time and condition, $\beta_2$ is equal to zero while $H_1$ states that $\beta_2$ is larger than zero.

$$H_0: \beta_2=0$$ 
$$H_1: \beta_2>0$$

## Description of files

- The file "Simulation.RMD" contains code for plotting results of the simulation study. This requires prior execution of the two other main functions "fct_data_generation_vec" and "fct_BayesianSSD".
- The file "fct_data_generation_vec" contains a vectorized version of the function to generate datasets under both hypotheses at hand and to calculate Bayes Factors for each dataset.
- The file "fct_BayesianSSD" contains the function needed to execute the sample size determination using the data generating function. 

## Algorithm 2
The file "Algorithm2" contains a refinement of the algorithm used in the simulation. Using a binary search, it reduces the amount of necessary iterations to maximally 12 (see Fu, Hoijtink, and Moerbeek, 2020 for a brief overview). 




