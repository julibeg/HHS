
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
import sys
sys.path.append(
    '/home/sogga/Dropbox/Julian_MSc_Project_Folder/Scripts/submission_code')
from final import *
# %% ##########################################################################
dists = squareform(np.arange(1, 11))
dists
# %% ##########################################################################
avg_pairwise_dist(dists, [0, 2, 3, 4])
# %% ##########################################################################
dists = pd.read_csv('test100.out', header=None)
dists.head()
dists.shape
# %% ##########################################################################
indices = np.array([0, 2, 3, 4, 66, 89, 95])
avg_pairwise_dist(dists.values, indices)
# %% ##########################################################################
indices = np.random.randint(0, 100, 54)
indices.sort()
indices = np.unique(indices)
avg_pairwise_dist(dists.values, indices)
indices
indices.shape
# %% ##########################################################################
gt = np.genfromtxt('200samples_10000snps.gt', delimiter=1)
phen = np.genfromtxt('200samples.phen', delimiter=1)
# %% ##########################################################################
gt[gt == 2] = np.nan
# %% ##########################################################################

PG_counts(gt.T[4762], phen)
# %% ##########################################################################
gt2 = np.array([np.nan, 1, 0, 0, 1, np.nan, 0])
phen2 = np.array([1, 0, 0, 1, 1, 0, 0])

PG_counts(gt2, phen2)
