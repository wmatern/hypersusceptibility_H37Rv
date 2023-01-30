#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import sklearn.metrics as skm

#Globals for plot settings
n_f = 0 #Global for setting figure numbers

def main():
    global n_f

    out_fold = 'output/'
    padj_thresh = 0.05
    lfc_thresh = 0.5
    label_thresh = 100.0
    save_plots = False
    run_plots = True

    #INH
    fname1 = out_fold+'summary_7d_INH.csv'
    fname2 = out_fold+'summary_14d_INH.csv'
    dose_col_names = ['H+++']
    union = compare_12h_vs_48h(fname1,fname2,padj_thresh,label_thresh,lfc_thresh,dose_col_names,run_plots,mode='union')
    if save_plots:
        plt.savefig(out_fold+'union_'+dose_col_names[-1]+'_7d_vs_14d'+'.svg')
    else:
        plt.show()

    inter = compare_12h_vs_48h(fname1,fname2,padj_thresh,label_thresh,lfc_thresh,dose_col_names,run_plots,mode='inter')
    inter.to_csv(out_fold+'intersection_INH.csv',index=False)
    #plot_12h_vs_48h(inh_sig, dose_col_names, label_thresh, run_plots)
    if save_plots:
        plt.savefig(out_fold+'inter_'+dose_col_names[-1]+'_7d_vs_14d'+'.svg')
    else:
        plt.show()

    #histogram
    if False:
        label = 'H+++'
        time = '_7d'
        plot_histograms(inh_merge,padj_thresh,label,time)

    #RMP
    fname1 = out_fold+'summary_7d_RMP.csv'
    fname2 = out_fold+'summary_14d_RMP.csv'
    dose_col_names = ['R+','R++','R+++']
    union = compare_12h_vs_48h(fname1,fname2,padj_thresh,label_thresh,lfc_thresh,dose_col_names,run_plots,mode='union')
    if save_plots:
        plt.savefig(out_fold+'union_'+dose_col_names[-1]+'_7d_vs_14d'+'.svg')
    else:
        plt.show()

    inter = compare_12h_vs_48h(fname1,fname2,padj_thresh,label_thresh,lfc_thresh,dose_col_names,run_plots,mode='inter')
    inter.to_csv(out_fold+'intersection_RMP.csv',index=False)
    #plot_12h_vs_48h(rmp_sig, dose_col_names, label_thresh, run_plots)
    if save_plots:
        plt.savefig(out_fold+'inter_'+dose_col_names[-1]+'_7d_vs_14d'+'.svg')
    else:
        plt.show()

    #histogram
    if False:
        label = 'R+++'
        time = '_7d'
        plot_histograms(inh_merge,padj_thresh,label,time)

    exit()

def plot_12h_vs_48h(df, dose_col_names, label_thresh, run_plots=True):
    plot_effect_sizes(df,dose_col_names[-1] + '_7d',dose_col_names[-1] + '_14d',label_thresh,run_plots)

def plot_effect_sizes(df, xlabel, ylabel, label_thresh, run_plots=True):
    global n_f
    x = df[xlabel].tolist()
    y = df[ylabel].tolist()
    rho,_ = st.pearsonr(x,y)
    spr,_ = st.spearmanr(x,y)
    x_01 = [u > 0 for u in x]
    y_01 = [v > 0 for v in y]
    phi = skm.matthews_corrcoef(x_01, y_01)
    if run_plots:
        fig, ax = plt.subplots()
        labels = df['uid'].tolist()
        plt.scatter(x,y)
        for i, txt in enumerate(labels):
            if (abs(x[i]) > label_thresh) and (abs(y[i]) > label_thresh):
                ax.annotate(txt,(x[i],y[i]))
        plt.xlim(-10*1.1,10*1.1)
        plt.ylim(-10*1.1,10*1.1)
        plt.xlabel(xlabel + ' Log2(Fold-Change)')
        plt.ylabel(ylabel + ' Log2(Fold-Change)')
        plt.axhline(0, color='black')
        plt.axvline(0, color='black')
        ax.annotate('Pearson:'+str(round(rho,2))+'\nSpearman:'+str(round(spr,2))+'\nPhi:'+str(round(phi,2)),xy=(.75,.1),xycoords='axes fraction')
        n_f += 1


def compare_12h_vs_48h(fname12h,fname48h,padj_thresh,label_thresh,lfc_thresh,dose_col_names,run_plots=True,mode='union'):
    clar12h = pd.read_csv(fname12h)[['uid','product']+dose_col_names+['pval','pval-adj (BH)']]
    clar12h = clar12h.rename({'pval-adj (BH)':'padj'},axis='columns')
    clar48h = pd.read_csv(fname48h)[['uid']+dose_col_names+['pval','pval-adj (BH)']]
    clar48h = clar48h.rename({'pval-adj (BH)':'padj'},axis='columns')

    print(dose_col_names[-1])
    
    clar_merge = clar12h.merge(clar48h, how='inner',on='uid',suffixes=('_7d','_14d'))
    if mode == 'union':
        clar_simp = clar_merge[((clar_merge['padj_7d'] < padj_thresh) & (clar_merge[dose_col_names[-1]+'_7d'].abs() > lfc_thresh)) | ((clar_merge['padj_14d'] < padj_thresh) & (clar_merge[dose_col_names[-1]+'_14d'].abs() > lfc_thresh))]
    elif mode == 'inter':
        clar_simp = clar_merge[((clar_merge['padj_7d'] < padj_thresh) & (clar_merge[dose_col_names[-1]+'_7d'].abs() > lfc_thresh)) & ((clar_merge['padj_14d'] < padj_thresh) & (clar_merge[dose_col_names[-1]+'_14d'].abs() > lfc_thresh))]
    else:
        exit('Unrecognized mode:',mode)

    plot_12h_vs_48h(clar_simp, dose_col_names, label_thresh, run_plots)

    return clar_simp

def plot_histograms(df,padj_thresh,label,time):
    hist_val = df[(df['padj'+time] < padj_thresh) & (df[label+time] < 0.0)][label+time]
    plt.hist(hist_val,bins=20)
    plt.xlabel(label+time)
    plt.show()

    clar_sig = df[(df['padj_12h'] < padj_thresh) & (df['padj_48h'] < padj_thresh)]
    plt.hist(clar_sig[clar_sig[label+time] < 0.0][label+time],bins=20)
    plt.xlabel(label+time)
    plt.show()

#TODO: Remove Functions Below
def get_sig_genes(filename,padj_thresh,dose_col_names):
    clar12h = pd.read_csv(filename)[['uid']+dose_col_names+['pval','pval-adj (BH)']]
    clar12h = clar12h.rename({'pval-adj (BH)':'padj'},axis='columns')
    clar12h = clar12h.sort_values(by=dose_col_names[-1],axis='index')
    clar12h = clar12h[clar12h['padj'] < padj_thresh]

    return clar12h

def get_all_genes(filename,padj_thresh,dose_col_names):
    clar12h = pd.read_csv(filename)[['uid']+dose_col_names+['pval','pval-adj (BH)']]
    clar12h = clar12h.rename({'pval-adj (BH)':'padj'},axis='columns')

    return clar12h

if __name__=='__main__':
    main()
