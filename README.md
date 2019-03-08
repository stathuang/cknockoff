# Conditional-Knockoff
This repository contains source code for generating conditional knockoffs (paper on [arXiv](https://arxiv.org/abs/1903.02806)).
Rendered tutorials demonstrating the usage of the code are available at the [page](http://lucasjanson.fas.harvard.edu/code/ConditionalKnockoffs).

- Folder `src`: R code for generating conditional knockoffs for 3 models
  - `cknock_ldg.R`: low-dimensional Gaussian models
  - `cknock_ggm.R`: Gaussian graphical models
  - `cknock_dgm.R`: discrete graphical models (note: the graph-expanding for general graph is unfinished)
  - `util.R`: utility functions
  - `knockoff_measure`: computing importance statistics, adapted from the R package `knockoff`
  
The R markdown files demonstrates the usage of the code and reproduce the tutorials. 
  - gaussian.Rmd: low-dimensional Gaussian models
  - ggm.Rmd: Gaussian graphical models
  - dgm.Rmd: discrete graphical models

