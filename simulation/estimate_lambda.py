import ges
import sempler
import numpy as np
import pandas as pd
import csv

## define helper functions
def binary(x, k=None):
    x = np.asarray(x).reshape(-1)
    base = 2

    if k is None:
        k = np.maximum(np.floor(np.log2(x)).astype(int) + 1, 1)
        if isinstance(k, np.ndarray):
            k = np.max(k)
    else:
        kmax = np.maximum(np.floor(np.log2(np.max(x))).astype(int) + 1, 1)
        assert k >= kmax, f"Provided k={k} is too small; need at least {kmax} bits."

    powers = base ** np.arange(k - 1, -1, -1)
    divs = np.floor_divide(x[:, None], powers)
    shifted = np.hstack([np.zeros((len(x), 1), dtype=int), base * divs[:, :-1]])
    r = divs - shifted

    return r.astype(int)

def CPDAG(lambda_est):
    lambda_fin = lambda_est
    d = len(lambda_est)
    for i in range(d):
        for j in range(i+1,d):
            if lambda_est[i,j] == 1 and lambda_est[j,i] == 1:
                lambda_fin[i,j] = -1
                lambda_fin[j,i] = -1
    return lambda_fin

def count_accuracy(B_true, B_est):
    """
    Original Source: https://github.com/30bohdan/latent-dag
    Args:
        B_true (np.ndarray): [d, d] ground truth graph, {0, 1}
        B_est (np.ndarray): [d, d] estimate, {0, 1, -1}, -1 is undirected edge in CPDAG

    # SHD = undirected extra (skeleton) + undirected missing (skeleton) + reverse (directed graph)
    # unoriented_correct = # undirected edges in the cpdag that has a corresponding true edge in the true dag
    """
    d = len(B_true)
    assert (len(B_est) == d)
    undirected_extra = 0
    undirected_missing = 0
    reverse = 0
    unoriented_correct = 0
    for i in range(d):
        for j in range(i + 1, d):
            undir_true = (B_true[i][j] == 1 or B_true[j][i] == 1)
            undir_est = (B_est[i][j] == 1 or B_est[i][j] == -1 or B_est[j][i] == 1 or B_est[j][i] == -1)

            if undir_true and (not undir_est):
                undirected_missing += 1
            elif (not undir_true) and undir_est:
                undirected_extra += 1
            elif undir_true and undir_est:
                if B_est[i][j] == -1 or B_est[j][i] == -1:
                    # Undirected edge in est
                    unoriented_correct += 1
                elif B_true[i][j] != B_est[i][j]:
                    # Directed edge in est, but reversed
                    reverse += 1
    return undirected_extra + undirected_missing + reverse
    # return {"shd": undirected_extra + undirected_missing + reverse,
    #         "undirected_extra": undirected_extra,
    #        "undirected_missing": undirected_missing,
    #        "reverse": reverse,
    #        "unoriented_correct": unoriented_correct}

## main code
pi = pd.read_csv('pi_collider_5000.csv', header = None) # change filename accordingly

Nrep = len(pi)
K = 3
num_samples = 2000

## change the true graph accordingly
# linear
# lambda_true = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
# collider
lambda_true = [[0, 1, 0], [0, 0, 0], [0, 1, 0]]
# dependent
# lambda_true = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]

err_shd = np.zeros(Nrep)

for i in range(Nrep):
    prop = pi.iloc[i]

    # generate pseudo-amples
    one_hot_samples = np.random.multinomial(n=1, pvals=prop, size=num_samples)
    samples = np.argmax(one_hot_samples, axis=1)
    data = binary(samples, K)

    # run GES
    lambda_est, score = ges.fit_bic(data)
    lambda_est = CPDAG(lambda_est)
    
    err_shd[i] = count_accuracy(lambda_true, lambda_est)

# estimation accuracy of Lambda
print(np.mean(err_shd))
