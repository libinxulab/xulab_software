"""
    train_prediction_model/onfig.py
    Dylan H. Ross

    Configuration parameters for training the CCS prediction model
"""

# pRNG seed for deterministic results in stochastic steps (e.g., splitting training/test set data)
seed = 1234

# GridSearch number of jobs
gs_n_jobs = 10

# number of individual clusters to split the dataset into
n_clusters = [4, 5]

# SVR hyperparameters C and gamma, and whether to permute them for each cluster
C = [1e3, 1e5]
gamma = [1e-3, 1e-1]
hyperparam_permute = True

