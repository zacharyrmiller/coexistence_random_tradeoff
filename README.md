Code to reproduce analysis and figures from "Coexistence of many species under a random competition-colonization trade-off" by Zachary R. Miller, Maxime Clenet, Katja Della Libera, François Massol, and Stefano Allesina. 

Email zachary.miller@yale.edu wth any questions.

Required packages: tidyverse, RColorBrewer, ggpubr, deSolve, truncnorm

This repository includes:

cc_functions.R: core simulation functions for the competition-colonization trade-off model, including functions to find the coexisting set of species given an initial pool and to obtain equilibrium occupancies for all species (without numerical integration)

distribution_definitions.R: functions to sample from example probability distributions used throughout the study

paper_figures_code.R: code to reproduce all figures in the main text of the paper

supplement_figures_code.R: code to reproduce figures in the supplementary information
