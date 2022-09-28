
### Predictive inference for intractable models via ABC

This repository contains the code to reproduce the simulation experiments of the paper **_"On predictive inference for intractable models via approximate Bayesian computation"_**. A preprint of the paper is available at: https://arxiv.org/abs/2203.12495.

The codebase also includes some additional features and simulation experiments (e.g. `simplemarkov.inference`). These are however not used in the paper and are not carefully tested.

To run the simulation experiments, the variable `opt$root` in `simplemarkov1D.inference.R` and in both `*_setup.R`-files need to be modified to point to the folder that contains the code files. Similarly, `opt$save.loc` needs to point to a folder where the computed results are to be saved.

It is also necessary to compile the C-codes of the simulation models using:

```
R CMD SHLIB simplemarkov_simul.c
R CMD SHLIB MG1_simul.c
R CMD SHLIB LV_simul.c
```
(The code has been tested using 64-bit Kubuntu 20.04 LTS with R version 3.6.3 "Holding the Windsock". On a Windows platform it might be necessary to additionally modify the `dyn.load`-command in each `*_inference.R`-file so that the resulting `*.dll`-files are loaded instead of the `*.so`-files used in Linux/OSX environments.)

The following R libraries are needed to run (some parts of) the code:

* matrixStats
* latex2exp
* coda

These libraries can be installed using `install.packages()`.
