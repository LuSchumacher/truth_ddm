# Cognitive Processes Underlying the Truth Effect

This repository contains the data and code for reproducing all results reported in our paper "Cognitive Processes Underlying the Repetition-Based Truth Effect: A Diffusion Decision Model Study".

## Content

- [Model implementation](model/truth_ddm.stan): Bayesian hierarchical Diffusion Decision Model (DDM) implementation in Stan.
- [Behavioral analysis](src/2_glm.R): Generalized linear model based data analysis.
- [Cognitive model estimation](src/3_model_fitting.R): Estimation of the DDM.
- [Cognitive model evaluation](src/4_model_evaluation.R): Evaluation of the DDM parameter estimates.
- [Posterior re-simulation](src/5_post_resimulation.R): Evaluation of the DDMs absolute fit to observed data via posterior re-simulations.
- [Parameter recovery study](src/6_param_recovery.R): Evaluation of the DDMs identifiability via parameter recovery study.

## Support

This work was supported by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation; grant number GRK 2277 ”Statistical Modeling in Psychology”)

## License

MIT
