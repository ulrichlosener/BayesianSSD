# Bayesian Sample Size Determination for Multilevel Models with Linear or Log-linear Growth

## File description:
 - **BayeSSD.R** performs Bayesian SSD for MLM. This requires the functions getbf()" and `getpower()`.
 - **getbf.R** simulates data and calculates the Bayes Factor. This is needed for running `BayeSSD()`.
 - **getpower.R** determines the power for a given sample size. This requires the function `getbf()` and is needed for running "BayeSSD()". This function can also be used to determine the power for a fixed N.
 - **simulation.R** executes the simulations necessary to produce the tables 2 and 3 in the paper.
 - **Old_functions** is a folder containing old versions of the functions above and others.

## Project description:
A priori sample size determination (SSD) is essential in designing cost-efficient trials and in avoiding underpowered studies. In addition, reporting a solid justification for a certain sample size is required by most ethical committees and many funding agencies. Most often SSD is based on null hypothesis significance testing (NHST), an approach that has recently received severe criticism. As an alternative, Bayesian evaluation of informative hypotheses has been developed. Informative hypotheses reflect specific theoretical and/or empirical expectations using (in)equality constraints on model parameters. Bayes factors quantify the relative support in the data for informative hypotheses (including the null hypothesis) without suffering from some of the drawbacks of NHST. SSD for Bayesian hypothesis testing relies on simulations and has only been studied recently. Available software for this is limited to simple models such as ANOVA and t-test, in which observations are assumed to be independent from each other. However, this assumption is rendered untenable when employing a longitudinal design where observations are nested within individuals. In that case, a multilevel model should be used. This paper provides researchers with a tool to perform SSD for multilevel models with longitudinal data in a Bayesian framework. To this end, we apply our methodology to various realistic scenarios and introduce an open source R function that enables researchers to tailor the simulation to their trial at hand. 

The paper concerning this Project has been submitted for publication to Behavior Research Methods.
This work is by Ulrich Lösener, PhD candidate at Utrecht University, supervised by Dr. ir. Mirjam Moerbeek and Prof. Dr. Herbert Hoijtink. This project is funded by an open-competition NWO grant obtained by Mirjam Moerbeek.
