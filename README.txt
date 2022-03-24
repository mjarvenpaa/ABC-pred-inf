
Code used for the paper "On predictive inference for intractable models via approximate Bayesian computation", preprint: https://arxiv.org/abs/2203.12495 

The code also includes some additional features and simulation experiments (e.g. 'simplemarkov.inference'). 
These are however not used in the paper and not carefully tested. 

To run the code it is necessary to compile the C-codes of the simulation models using:

R CMD SHLIB MG1_simul.c
R CMD SHLIB LV_simul.c
R CMD SHLIB simplemarkov_simul.c

(The code has been tested using 64-bit Kubuntu 20.04 LTS with R version 3.6.3 "Holding the Windsock". 
On a Windows platform it might be necessary to modify the 'dyn.load'-commands so that the resulting *.dll-files are loaded instead of the *.so-files used in Linux/OSX environments.)

Furthermore, 'opt$root' variable e.g. in MG1_setup.R needs to be modified to point to the folder that contains the code files. 
Similarly 'opt$save.loc' needs to point to a folder where the computed results are to be saved. 

The following two R libraries are needed to run (some parts of) the code:
- matrixStats
- latex2exp
These can be installed in a usual way with 'install.packages'.
