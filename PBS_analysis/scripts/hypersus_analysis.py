#!/usr/bin/env python3
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from JT_test import JT_test_pd
from hypersus_helper import norm_data_for_JT, estimate_lfc, summary_data

#parameters
plots_on = True
run_calcs = True
out_fold = 'output/hypersus_analysis/'
pseudocount = 4

#Seed the random number generator - ensure consistent results each time
np.random.seed(0)

def main():
    #Initializations
    n_f = 0 #Figure index (adds one after each figure plot)

    #Read in and organized data, annot_df contains annotations of each site.
    full_df,full_annot_df = read_in_tn_data()

    ##Isolate single gene names in uid_df (for merging back after processing)
    uid_df = full_annot_df[['uid','product']].drop_duplicates(subset='uid') #Gene by gene annotations
    uid_df = uid_df[~uid_df['uid'].isna()]
    func = lambda x: x.find(';') == -1
    uid_df = uid_df[uid_df['uid'].map(func)]
    uid_df = uid_df.set_index('uid')

    #Remove sites with low read counts in the input
    lowrm_df = full_df[(full_df['0d']['InB'].sum(axis=1) > 10)] 
    annot_df = full_annot_df[(full_df['0d']['InB'].sum(axis=1) > 10)]

    #norm_7d_df = norm_df['7d']
    #norm_14d_df = norm_df['14d']

    lfc_7d_df = estimate_lfc(lowrm_df['7d'],annot_df,pseudocount)
    norm_7d_df = norm_data_for_JT(lowrm_df['7d'],annot_df)

    lfc_14d_df = estimate_lfc(lowrm_df['14d'],annot_df,pseudocount)
    norm_14d_df = norm_data_for_JT(lowrm_df['14d'],annot_df)

    #Plot histogram of log fold change
    #plt.hist(lfc_48h_df['C++'].dropna(),bins=50, color='r')
    #plt.show()

    if run_calcs:
        # Computations using 7d timepoint only
        treatment = ['ND','H+++']
        #merged_df = summary_data_R(norm_7d_df,lg2_mmfc_7d_df,treatment,N_perm)
        merged_df = summary_data(norm_7d_df,lfc_7d_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid').drop(['ND'],axis=1)
        merged_df.to_csv(out_fold+'summary_7d_INH.csv')

        treatment = ['ND','R+','R++','R+++']
        #merged_df = summary_data_R(norm_7d_df,lg2_mmfc_7d_df,treatment,N_perm)
        merged_df = summary_data(norm_7d_df,lfc_7d_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid').drop(['ND'],axis=1)
        merged_df.to_csv(out_fold+'summary_7d_RMP.csv')

        #Using 14d time point
        treatment = ['ND','H+++']
        #merged_df = summary_data_R(norm_14d_df,lg2_mmfc_14d_df,treatment,N_perm)
        merged_df = summary_data(norm_14d_df,lfc_14d_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid').drop(['ND'],axis=1)
        merged_df.to_csv(out_fold+'summary_14d_INH.csv')

        treatment = ['ND','R+','R++','R+++']
        #merged_df = summary_data_R(norm_14d_df,lg2_mmfc_14d_df,treatment,N_perm)
        merged_df = summary_data(norm_14d_df,lfc_14d_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid').drop(['ND'],axis=1)
        merged_df.to_csv(out_fold+'summary_14d_RMP.csv')

def read_in_tn_data():
    #Read in and organized data. Return a dataframe containing organized raw read count data.

    #Read in csv files of counts
    csv1 = 'output/fastq_back.csv'
    
    full_df = pd.read_csv(csv1,dtype={'regulatory_class':str,'bound_moiety':str})
    annot_df = full_df.iloc[:,0:6].rename(columns={'unique_identifier (locus_tag or record_id_start_end_strand)':'uid'})

    mindex = pd.MultiIndex.from_arrays([annot_df['contig'],annot_df['insertion_site']])
    full_df = full_df.drop(labels=['contig','insertion_site','unique_identifier (locus_tag or record_id_start_end_strand)','product','regulatory_class','bound_moiety'],axis=1)
    full_df.index = mindex
    annot_df=annot_df.iloc[:,2:]
    annot_df.index = mindex
    #full_df.to_csv(out_fold+'hypersus_testing.csv')

    labels_df = pd.read_csv("input/Sample_Labels.csv")
    labels_df = labels_df[['Sample Description','SRA Accession Number']]
    labels_df = labels_df.rename(columns={'Sample Description':'ID','SRA Accession Number':'PrepLabel'})

    #Read in data and rename for ease of reference
    samp_column_names = list()
    tuples = list()
    for i in range(len(labels_df["ID"])):
        old_name = 'read_count (' + str(labels_df["PrepLabel"][i]) + ')'
        if labels_df["ID"][i][0:3] == 'InA':
            tuples.append(('-17d','InA'))
        elif labels_df["ID"][i][0:3] == 'InB':
            tuples.append(('0d','InB'))
        elif (labels_df["ID"][i][0:2] == 'd7'):
            tuples.append(('7d',labels_df["ID"][i][3:]))
        elif (labels_df["ID"][i][0:3] == 'd14'):
            tuples.append(('14d',labels_df["ID"][i][4:]))
        else:
            print(labels_df["ID"][i])
            exit("Unexpected ID")
        samp_column_names.append(old_name)
    
    #Build a multi-Index for the columns - You want to be able to refer to each group: (d7,d14) and (H+,R+++...)
    full_df = full_df[samp_column_names]
    mindex = pd.MultiIndex.from_tuples(tuples, names=['time','abx'])
    full_df.columns=mindex

    #Sort the column names for convenience
    full_df = full_df.sort_index(axis=1,level=['time','abx'])
    #pd.concat([annot_df,full_df],join='inner',axis=1,levels=['time','abx']).to_csv(out_fold+'Organized_raw_data.csv')

    return(full_df, annot_df)

if __name__ == '__main__':
    main()
