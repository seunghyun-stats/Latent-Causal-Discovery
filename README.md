# Latent-Causal-Discovery

This repository contains the MATLAB and Python codes used for experiments in the manuscript [Lee & Gu (2025), Identifiability of latent causal graphical models without pure children]. The codes estimate the observed-to-latent bipartite graph and the latent causal graph. The main algorithm in this code is based on Algorithm 1 in the paper by Ma, Ouyang & Xu (2023) (https://link.springer.com/article/10.1007/s11336-022-09867-5).

### For simulations:
To run simulations, go to the folder `simulations` and run the `simulation_main.m` file. The simulation settings such as the true graphical structures and conditional probabilities can be changed within the script. To estimate the latent DAG $\Lambda$, run the `estimate_lambda.py` file (in Python).

### Additional codes:
The folder `utilities` contains helper functions required to implement the main algorithms. We recommend adding this folder using the `add path 'utilities'` command in MATLAB.

### Dependencies:
MATLAB
- 9.13.0.2126072 (R2022b)

Python
- numpy 1.26.2
- pandas 2.2.2
- ges 1.1.1
