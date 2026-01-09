# STICCC Analysis Scripts

These are the R scripts used to produce the analysis in the manuscript "Dissecting reversible and irreversible single cell state transitions from gene regulatory networks"


The STICCC software package is required to run these scripts and is available here: https://github.com/lusystemsbio/viccc/

The scripts contain the following content:

* ..._Inference.R: simulations of synthetic circuits and subsequent application of STICCC
* ..._trajectory_validation.R: comparisons between STICCC predictions and simulated trajectories
* ...RNAVelocity_Comparison.R: comparisons between STICCC predictions and RNA velocity predictions
* CTS_dropout_comparison.R: effect of measurement noise on STICCC predictions
* CTS_signaling.R: STICCC applied to multiple timepoints of a simulated dataset under perturbation
* Cellcycle_U2OS_comparison.R: additional application to a cell cycle dataset
* Repressilator_downsample.R: effect of dataset size on STICCC predictions
* sticcc_analysis_utilities.R: helper functions for some analysis and visualization of STICCC results beyond what the package includes

For additional information, contact the authors or see the original paper: https://www.biorxiv.org/content/10.1101/2024.08.30.610498v1

