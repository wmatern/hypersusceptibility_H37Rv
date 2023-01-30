import pandas as pd
import numpy as np
import statsmodels.stats.multitest as smm
import time
from scipy import stats
from JT_test import JT_test_pd

def estimate_lfc(count_df,annot_df,pseudocount,ctrl_grp='ND'):
    # Compute approximate lfc after adding pseudocount
    df = count_df + pseudocount
    norm_df = df / df.sum(axis=0)

    #Remove sites without features
    norm_df = norm_df[~annot_df['uid'].isna()]
    annot_df = annot_df[~annot_df['uid'].isna()]

    #Remove sites which contain more than one feature.
    func = lambda x: x.find(';') == -1
    norm_df = norm_df[annot_df['uid'].map(func)]
    annot_df = annot_df[annot_df['uid'].map(func)]

    #Compute Effect sizes (LFC)
    mnorm_df = norm_df.groupby(level=0,axis=1,sort=False).mean()
    lg2_mfc_df = mnorm_df.div(mnorm_df[ctrl_grp],axis=0).apply(np.log2)#Mean log fold change

    #Take the median of the log fold change across TA sites in the same gene
    lfc_df = lg2_mfc_df.groupby(by=annot_df['uid'],sort=False).median()
    
    return lfc_df

def norm_data_for_JT(count_df,annot_df):
    # Normalize and clean data in preparation to pass to JT
    norm_df = count_df / count_df.sum(axis=0)

    #Remove sites without features
    norm_df = norm_df[~annot_df['uid'].isna()]
    annot_df = annot_df[~annot_df['uid'].isna()]

    #Remove sites which contain more than one feature.
    func = lambda x: x.find(';') == -1
    norm_df = norm_df[annot_df['uid'].map(func)]
    annot_df = annot_df[annot_df['uid'].map(func)]

    norm_df = annot_df[['uid']].merge(norm_df,how='inner',left_index=True,right_index=True)
    norm_df = norm_df.set_index('uid')

    return norm_df

def summary_data(norm_12h_rm0_df,mmfc_12h_df,treatment):
    #Runs JT testing, multiple hypothesis correction, and returns summarized data
    t0 = time.time()
    df = norm_12h_rm0_df[treatment]
    df_index = df.index
    df = df.reset_index(drop=True)

    df['pval'] = JT_test_pd(df,treatment,direction='lessthan',random_smooth=True)
    df.index = df_index

    p_ls = []
    name_ls = []
    for name,g in df.groupby(['uid'],sort=False):
        name_ls.append(name)
        Zscore,_ = stats.combine_pvalues(g['pval'],method='stouffer')
        p_pooled = stats.norm.sf(abs(Zscore))*2
        p_ls.append(p_pooled)
    _,p_adj_ls,_,_ = smm.multipletests(p_ls, alpha=0.05, method='fdr_bh')
    pval = pd.DataFrame({'uid':name_ls,'pval':p_ls,'pval-adj (BH)':p_adj_ls})

    merged_df = mmfc_12h_df[treatment].merge(pval,how='left',left_on='uid',right_on='uid')
    merged_df.index = merged_df['uid']
    merged_df = merged_df.drop('uid',axis=1)
    print('JT test run time:', time.time()-t0)
    return merged_df
