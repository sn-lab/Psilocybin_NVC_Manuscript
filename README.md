# Psilocybin_NVC_Manuscript

Repository for documentation, code, and smaller datasets associated with the SN-lab manuscript on psilocybin-induced changes in mouse NVC. Each section below details the contents of a folder containing processed data tables and/or code for analyzing and plotting the data within those tables.

## Two-photon (2P) Imaging Analysis

To be added.

## Widefield Imaging Analysis

To be added.

## BOLD Simulations

The code used to simulate single-node and whole-brain neural mass model activity, blood flow, and BOLD were modifed based on the [neurolib repository](https://neurolib-dev.github.io/). After installing neurolib based on the instructions in that library, four files were manually modified to enable custom adjustments of Balloon-Windkessel NVC model parameters. The files which were modified can be found in the [neurolib/neurolib/models](https://github.com/neurolib-dev/neurolib/tree/master/neurolib/models) folder of the original repository and the "aln" and "bold" subfolders within that directory. The modified files are contained in this "Bold Simulations folder, in the "custom files" subfolder.

**To perform whole-brain simulations:** After the neurlib library was installed and the relevent files were replaced with the modified files, simulations were performed using the "bold_sim.py" script in python. This script contains three code sections for running BOLD simulations. In the "run_example" section, a single whole-brain simulation is performed, and the resulting neural activity data, BOLD data, and functional connectivity (FC) matrix is saved. In the run_series section, a set of simulations is performed under identical NVC parameters but with different random seeds to initiate the model with different neural activity patterns, and the resulting FC matrix is saved for each. In the "run_matrix" section, a set of simulations are repeated across a matrix of "Gamma" and "K" NVC parameters values, and the global FC strength of each simulation is saved. To plot the resulting data from these simulations, the script "analyze_bold_80nodes.m" was run in Matlab.

**To perform single-node simulations:** The Balloon-Windkessel model from neurolib was translated for Matlab for both simulation and parameter fitting with a square-wave neural activity input. The "integrateBOLD" function was used to simulate the flow and BOLD signals from a square wave activity input. The "fitBOLD" function was used to fit measured vessel flow responses (stored in the "vel_avgs_2p.mat" and "speed_avgs_widefield.mat" files) to a simulated flow response from optimized NVC parameters. All of these single-node simulation functions were called, and all plots were generated, from the "analyze_bold_1node.m" script.
