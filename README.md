# Latent-Causal-Discovery

This repository contains the MATLAB and Python codes for the manuscript [Lee, S. and Gu, Y. (2025), Identifiability of latent causal graphical models without pure children]. The main algorithm is based on Algorithm 1 in the paper by Ma, Ouyang & Xu (2023) (https://link.springer.com/article/10.1007/s11336-022-09867-5).

### For simulations:
To run simulations, go to the folder `simulations` and run the `simulation_main.m` file. The simulation settings such as the true graphical structures and conditional probabilities can be changed within the file. To estimate the latent DAG $\Lambda$, run the `estimate_lambda.py` file.

### Additional codes:
The folder `utilities` contains helper functions that will be required to implement the main algorithms. We recommend adding this folder using the `addpath 'utilities'` command in MATLAB.

