from julsy.stuff import parallel_apply_along_axis, na_bincount,\
    zero_D_arrs_to_df
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
import sys
sys.path.append(
    '/home/sogga/Dropbox/Julian_MSc_Project_Folder/Scripts/submission_code')
from final import *
from adapt_cross_resist import adapt_iter, adapt_cross_resistance
from tqdm import tqdm
# %% ##########################################################################
gt = np.genfromtxt('test.gt', delimiter=1)
gt[gt == 2] = np.nan
phen = np.genfromtxt('test.phen', delimiter=1)
dists = pd.read_csv('test.dists', header=None)
# %% ##########################################################################
gt = np.genfromtxt('9586snps_199samples.gt.filt', delimiter=1).T
gt[gt == 2] = np.nan
phen = np.genfromtxt('199samples.phen.filt', delimiter=1)
dists = pd.read_csv('199samples_9586snps.dists.filt', header=None)
# %% ##########################################################################
phen = pd.Series(phen)
gt = pd.DataFrame(gt)
# %% ##########################################################################
print(dists.shape)
print(phen.shape)
print(gt.shape)
# %% ##########################################################################


def full(gt, phen, avg_pw_dists, gt_weights='SNP', rel_gt_weight=None,
         phen_weights=None, rel_phen_weight=None, P0G1_extra_weight=1,
         dist_weight=1, df=None, P1G1_filter=2, **adapt_kwargs):
    """
    full HHS function; requires at least a GT matrix, phenotype vector and
    average pairwise distances.
    """

    if isinstance(phen, pd.DataFrame):
        phen = phen[phen.columns[0]]
    p1g1, p0g0, p1g0, p0g1 = PG_counts(gt, phen)
    if phen_weights is None:
        phen_weights = phen.shape[0] / 2 / np.bincount(phen)
    if gt_weights == 'SNP':
        gt_weights = gt.apply(
            lambda col: (col.shape[0] - col.isna().sum()) / 2 /
            col.value_counts()).values
    elif gt_weights == 'all':
        gt_weights = (np.prod(gt.shape) - gt.isna().sum().sum()) / \
            2 / na_bincount(gt.values.flatten())
    if rel_gt_weight is not None:
        gt_weights = rel_gt_weight * gt_weights + (1 - rel_gt_weight)
    if rel_phen_weight is not None:
        phen_weights = rel_phen_weight * phen_weights + (1 - rel_phen_weight)
    counts = p1g1 * phen_weights[1] * gt_weights[1] - \
        p0g1 * phen_weights[0] * gt_weights[1] * P0G1_extra_weight
    scores = counts + (avg_pw_dists - 1) * counts * dist_weight

    scores[p1g1 < P1G1_filter] = 0
    scores[scores < 0] = 0

    print(phen_weights)
    print('---------')
    print(gt_weights.sum())
    print(gt_weights.max())
    print('---------')
    print(p1g1.sum())
    print(p1g1.max())
    print('---------')
    print(p0g1.sum())
    print(p0g1.max())
    print('---------')
    print(counts.sum())
    print(counts.max())
    print('---------')
    print(scores.sum())
    print(scores.max())
    print('---------')

    scores_adapted = adapt_iter(
        gt, phen, scores, include_orig_scores=True, **adapt_kwargs).T

    if df is None:
        df = pd.DataFrame(index=gt.columns)
    df_2 = zero_D_arrs_to_df(arrs=(p1g1, p0g1, p1g0, p0g0, counts),
                             columns=('p1g1', 'p0g1', 'p1g0', 'p0g0', 'counts'),
                             index=df.index).astype(
        {x: int for x in ('p1g1', 'p0g1', 'p1g0', 'p0g0')})
    scores_adapted.columns = [f'scores{x}' for x in scores_adapted.columns]
    return pd.concat((df, df_2, scores_adapted), axis=1)


def adapt_iter(gt_orig, phen, scores, n_iterations=10, effect_size=0.001,
               final_rescale=True, save_every=100, include_orig_scores=True,
               P0G1_weight=0, threshold=1e-10, intermittent_rescale=False):
    """
    iterative HHS elimination
    """
    print(scores.max())
    orig_scores = scores.copy()
    orig_scores[orig_scores < threshold] = 0
    gt = gt_orig.copy()
    gt[np.isnan(gt)] = 0
    saved_iterations = np.arange(0, n_iterations, save_every) + 1
    if saved_iterations[-1] != n_iterations:
        saved_iterations = np.append(saved_iterations, n_iterations)
    score_arr = np.empty((len(saved_iterations), len(scores)))
    for i in tqdm(range(int(n_iterations))):
        scores = adapt_cross_resistance(gt, phen, scores, effect_size,
                                        rescale=intermittent_rescale,
                                        P0G1_weight=P0G1_weight,
                                        threshold=threshold)
        if i % save_every == 0:
            saved_idx = i // save_every
            score_arr[saved_idx] = scores
    score_arr[-1] = scores
    print(scores.max())
    if final_rescale:
        score_arr *= orig_scores.sum() / score_arr.sum(1).reshape(-1, 1)
    score_df = pd.DataFrame(score_arr, index=saved_iterations,
                            columns=gt.columns)
    if include_orig_scores:
        score_df = pd.concat((pd.DataFrame(orig_scores.reshape(1, -1),
                                           index=[0],
                                           columns=score_df.columns),
                              score_df), axis=0)
    return score_df

# %% ##########################################################################


# avg_dists = all_relative_avg_dists(gt, phen, dists)
# avg_dists.shape

# %% ##########################################################################
%time res = full(gt, phen, avg_dists, gt_weights='all', rel_phen_weight=0.2, rel_gt_weight=0.2, n_iterations=1e4, effect_size=1e-4, save_every=10)
# %% ##########################################################################
last = res.iloc[:, -1]
last.loc[last > 0]
# %% ##########################################################################
res.sum().iloc[5:]
# %% ##########################################################################
(res > 0).sum().iloc[5:]
# %% ##########################################################################
res.max().iloc[5:]
