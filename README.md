# Stochastic SMAD simulation code
Simulation code of the TGF-β/SMAD pathway in single cells as used for the paper
> Data-based stochastic modelling reveals sources of activity bursts in single-cell TGF-β signaling

by N Kolbe, L Hexemer, LM Bammert, A Loewer, M Lukáčová-Medviďová and S Legewie

A preprint of the paper is available on [arXiv](https://arxiv.org/abs/2107.11770).

Matlab codes were developed by N Kolbe, L Hexemer, LM Bammert and data contained in this Repository
was produced by A Loewer and J Strasen and partly previously published on [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.hc5dp).

## Contents
The root folder contains two scripts with included documentation to simulate the stochastic SMAD model 1) for a single cell (`single_path.m`) and 2) for a population of cells (`many_paths.m`). The second script allows for a comparison to experimental data using the objective function defined in the paper. In addition, the script `tune_burst_detection.m` shows how the parameters of the burst detection were optimized in the paper.

The functions in the folder `+burstDetection` perform single cell analysis and compute the objective function. The folder `+forward` contains the SDE solver for the model and its dependencies. Experimental data is located in the folder `data`. The contents of all three folders are needed to run the scripts.

## Usage
Copy or clone the full content of this repository to your computer. The code can then be run in MATLAB from within the main folder of the copied code. The code is commented and you can get easily started by going through the scripts. 

The code was developed in MATLAB 2017a and it will likely be working for newer versions (it was verified to work in MATLAB 2020b). The programs make use of MATLAB packages and for the burst analysis the Signal Processing Toolbox is needed. 

For the parameter estimation in the paper, this code was used together with the global optimization toolbox [MEIGO](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-136) available [here](https://bitbucket.org/jrbanga_/meigo64/src/master/).

In case of questions or problems please contact the authors of the paper or file a GitHub issue for this repository.
