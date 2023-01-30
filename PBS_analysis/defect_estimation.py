#!/usr/bin/env python3
import hypersus_analysis as ha
import pandas as pd
import matplotlib.pyplot as plt

#
#import pandas as pd 
#import matplotlib.pyplot as plt
#import math
#import statsmodels.stats.multitest as smm
#import numpy as np
#from scipy import stats
#import time
#
##packages for running r code
#from rpy2.robjects import pandas2ri
#pandas2ri.activate()
#from rpy2.robjects.packages import importr
##Import R gMWT package for fast JT computation
#gmwt = importr("gMWT")
#
##these are constants used for plotting on log scale
#MAXLOG = 10
#MINLOG = -10

#parameters
plots_on = True
run_calcs = True
out_fold = 'output/'

def main():
    #Initializations
    n_f = 0 #Figure index (adds one after each figure plot)

    #Read in and organized data, annot_df contains annotations of each site.
    full_df,annot_df = ha.read_in_tn_data()

    #Normalize to sum of the reads in each sample. #TODO: For samples with loss of diversity you should use a different normalizer - some estimate of the wildtype. 
    norm_df = full_df / full_df.sum(axis=0)

    #Remove sites which have no reads in Input samples.
    sum_temp_sr = (full_df['-14d']['InA'].sum(axis=1) > 0)
    norm_df = norm_df[sum_temp_sr]
    annot_df = annot_df[sum_temp_sr]

    #Change the index of the dataframes to use the uid (gene names)
    norm_df,annot_df = ha.reindex_by_uid(norm_df,annot_df)
    
    norm_dmso_df = pd.concat([norm_df['0d']['InB'],norm_df['7d']['ND'],norm_df['14d']['ND']],axis=1,sort=False)
    norm_dmso_df.columns = ['In']*3 + ['12h']*3 + ['48h']*3

    #Calculate mean/median FC (DMSO12h,DMSO48h:In) for each site.
    mnorm_df = norm_dmso_df['In'].mean(axis=1)
    mfc_df = norm_dmso_df.groupby(level=0,axis=1,sort=False).mean().div(mnorm_df,axis=0) #Mean fold change
    mmfc_df = mfc_df.groupby(level=0,axis=0,sort=False).median()

    #Log scaled fold change
    lg2_mmfc_df = mmfc_df.applymap(ha.lg2plus)

    #Statistical Computation
    if run_calcs:
        treatment = ['In','12h','48h']
        merged_df = ha.summary_data_R(norm_dmso_df,lg2_mmfc_df,treatment,3)
        merged_df.to_csv(out_fold+'summary_DMSOvsIn.csv')

    if plots_on:
        namex = '12h'
        namey = '48h'
        plt.figure(n_f)
        plt.scatter(lg2_mmfc_df[namex],lg2_mmfc_df[namey])
        plt.xlim(xmin=ha.MINLOG*1.1,xmax=ha.MAXLOG*1.1)
        plt.ylim(ymin=ha.MINLOG*1.1,ymax=ha.MAXLOG*1.1)
        plt.xlabel('Log2-RelativeFitness_'+namex)
        plt.ylabel('Log2-RelativeFitness_'+namey)
        #plt.show()
        plt.savefig(out_fold+'Fitness_'+namex+'_vs_'+namey+'.png')
        n_f+=1

if __name__ == '__main__':
    main()
