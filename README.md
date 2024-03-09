# Bayesian Causal Synthesis for Meta-inference on Heterogeneous Treatment Effects 

This repository provides R code implementing Bayesian causal synthesis (BCS) for combining multiple estimators of heterogeneous treatment effects (HTE), as proposed by the following paper.

Sugasawa, S., Takanashi, K., McAlinn, K. and Airoldi, E. M. (2023). Bayesian Causal Synthesis for Meta-Inference on Heterogeneous Treatment Effects. ([arXiv:2304.07726](https://arxiv.org/abs/2304.07726))

The repository includes the following files.

- `BCS.R`: R code implementing BCS given outputs of multiple HTE estimators
- `BCS-default.R`: R code implementing default version of BCS combining three estimators of HTE (Bayesian causal forest, linear models and additive models) 
- `demo-ACIC.R`: Example using simulated dataset of Atlantic causal inference conference data analysis challenge 2017 (ACIC)
- `demo-Sim.R`: Example using simulated dataset under RCT 
- `Covariate.RData`: Covariate data used in `demo-ACIC.R`

Please contact ``sugasawa (at) econ.keio.ac.jp`` if you have any questions. 
