<<<<<<< HEAD

# Bayesian Sample Size Determination for Multilevel Models with Longitudinal Data

The file "Simulation.RMD" contains simulated multilevel data to which different multilevel models are fitted and estimated using Bayesian and frequentist estimation techniques.

## Data generation

In the trajectory of treatment interventions individuals may vary both in their initial symptom level and rate of change over time. That is, people might start out the trial with different severity of their symptoms and some of them might get better while others stay the same or even get worse. Because of random allocation to treatment conditions, the initial (pre-treatment) differences in symptom levels is unrelated to the treatment groups. The further symptom trajectory of an individual, however, depends on which treatment condition they are assigned to. Therefore, I choose a model which accounts for interindividual differences at the first measurement (random intercept) that cannot be explained by the treatment condition and which allows individuals to have different rates of change over time (random slopes). These differences in slopes, in turn, can be (partly) explained by membership of a treatment condition. 

It is common in longitudinal intervention data that measurements closer in time correlate higher with each other than measurements further separated in time. Also, variability in the outcome typically increases over time as individuals tend to be more similar in the beginning of a study and subsequently change at different rates. This renders the assumption of compound symmetry (constant variances and covariances over time) practically untenable. To account for this, a multilevel regression model with random intercept and random slopes is used in the following simulation.  

Level 1 (within subjects) regression equation: 
$$y_{ij} = b_{0i} + b_{1i}t_j + e_{ij}$$ where $y_{ij}$ denotes the score on the outcome variable of individual i ($i=1,2,...,N$) at measurement occasion j ($j=1,2,...,n$) and $$e_{ij} \sim N(0, \sigma^2)$$

Level 2 (between subjects) regression equations: 
$$b_{0i} = \beta_0 + u_{0i}$$ 
$$b_{1i} = \beta_1 + \beta_2 X_i + u_{1i}$$ where $X_i$ denotes the treatment condition, $u_{0i}$ is the individual deviation from the treatment group intercept, and $u_{1i}$ is the individual deviation from the treatment group slope. Note that membership of the treatment group ($X_i$) does not predict initial symptom levels ($b_{0i}$) but only the rate of change over time ($b_{1i}$).

By substitution we get the total regression equation containing both levels.

$$y_{ij} = \beta_0 + u_{0i} + (\beta_1 + \beta_2 X_i + u_{1i})*t_j + e_{ij}$$
Or
$$y_{ij} = \beta_0 + u_{0i} + \beta_1 t_j + \beta_2 X_i t_j + u_{1i} t_j + e_{ij}$$

Here, $\beta_0$ represents the average value for y at $t_0$ and $\beta_1$ the average rate of change for subjects in the $X_0$ condition. The average difference in the rate of change of y scores of subjects in the $X_0$ condition (relative to those in the $X_1$ condition) is represented by $\beta_2$. Thus, $\beta_2$ denotes the interaction between time and treatment, indicating the difference in slopes between the two treatment groups. To test whether a treatment intervention has an effect on symptom level, this is the parameter of interest. 

The population distribution of $u_{0i}$ and $u_{1i}$ is assumed to be bivariate normal $N(0, \Sigma_u)$ with 
$$\Sigma_u = \begin{bmatrix} \sigma^2_{u0} & \sigma_{u0u1}\\
\sigma_{u0u1} & \sigma^2_{u1}
\end{bmatrix}$$


969a27a667c9dd6ef78e19807e057d7212aa4bc2
