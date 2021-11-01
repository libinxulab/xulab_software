"""
    train_prediction_model/onfig.py
    Dylan H. Ross

    Configuration parameters for training the CCS prediction model
"""

# pRNG seed for deterministic results in stochastic steps (e.g., splitting training/test set data)
seed = 1234

# number of individual clusters to split the dataset into
n_clusters = [5]

# SVR hyperparameters C and gamma
C = [100, 1000]
gamma = [0.001, 0.01]
