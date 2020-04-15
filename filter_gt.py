from julsy.stuff import parallel_apply_along_axis, na_bincount,\
    zero_D_arrs_to_df
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
import sys
sys.path.append(
    '/home/sogga/Dropbox/Julian_MSc_Project_Folder/Scripts/submission_code')
from final import *
# %% ##########################################################################
gt = np.genfromtxt('200samples_10000snps.gt', delimiter=1)
gt[gt == 2] = np.nan
phen = np.genfromtxt('200samples.phen', delimiter=1)
dists = pd.read_csv('200samples_10000snps.dists', header=None)
# %% ##########################################################################
gt = pd.DataFrame(gt)
phen = pd.Series(phen)
# %% ##########################################################################
gt_filt = filter_GT_arr(gt, SNP_NA_threshold=0.05, strain_NA_threshold=0.05)
phen_filt = phen.loc[gt_filt.index]
dists_filt = dists.loc[gt_filt.index, gt_filt.index]
dists_filt.shape
gt_filt.shape
# %% ##########################################################################
gt_filt.fillna(2, inplace=True)
np.savetxt('9586snps_199samples.gt.filt', gt_filt.values.T.astype(int),
           delimiter='', fmt='%g')
np.savetxt('199samples.phen.filt', phen_filt.values.astype(int),
           delimiter='', fmt='%g', newline='')
dists_filt.to_csv('199samples_9586snps.dists.filt', index=False, header=False)
