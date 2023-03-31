#!/usr/bin/env python3
import pandas as pd 
import numpy as np
from hypersus_helper import norm_data_for_JT, estimate_lfc, summary_data

#parameters
run_calcs = True
out_fold = 'output/hypersus_analysis/'
pseudocount = 4

#Seed the random number generator - ensure consistent results each time
np.random.seed(0)

def main():
    #Initializations
    n_f = 0 #Figure index (adds one after each figure plot)

    #Read in and organize data, annot_df contains annotations of each site.
    full_df,full_annot_df = read_in_tn_data()

    ##Isolate single gene names in uid_df (for merging back after processing)
    uid_df = full_annot_df[['uid','product']].drop_duplicates(subset='uid') #Gene by gene annotations
    uid_df = uid_df[~uid_df['uid'].isna()]
    func = lambda x: x.find(';') == -1
    uid_df = uid_df[uid_df['uid'].map(func)]
    uid_df = uid_df.set_index('uid')

    lowrm_df = full_df[(full_df['0d']['In'].sum(axis=1) > 10)] #Remove sites with low read counts in the input
    annot_df = full_annot_df[(full_df['0d']['In'].sum(axis=1) > 10)]

    lfc_6d_df = estimate_lfc(lowrm_df['6d'],annot_df,pseudocount)
    norm_6d_df = norm_data_for_JT(lowrm_df['6d'],annot_df)

    if run_calcs:
        # Computations using 6d timepoint
        treatment = ['ND','H+']
        #merged_df = summary_data_R(norm_6d_df,lg2_mmfc_6d_df,treatment,N_perm)
        merged_df = summary_data(norm_6d_df,lfc_6d_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df = merged_df.drop(['ND'],axis=1)
        merged_df.to_csv(out_fold+'summary_6d_INH.csv')

        treatment = ['ND','R+']
        #merged_df = summary_data_R(norm_6d_df,lg2_mmfc_6d_df,treatment,N_perm)
        merged_df = summary_data(norm_6d_df,lfc_6d_df,treatment)
        merged_df = uid_df.merge(merged_df,how='left',on='uid')
        merged_df = merged_df.drop(['ND'],axis=1)
        merged_df.to_csv(out_fold+'summary_6d_RMP.csv')

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

    #labels_df = pd.read_excel("tn-seq_input/TnPrep_Labels_7H9.xls",nrows=18)
    labels_df = pd.read_csv("input/Sample_Labels.csv")
    labels_df = labels_df[['Sample Description','SRA Accession Number']]
    labels_df = labels_df.rename(columns={'Sample Description':'ID','SRA Accession Number':'PrepLabel'})

    #Read in data and rename for ease of reference
    samp_column_names = list()
    tuples = list()
    for i in range(len(labels_df["ID"])):
        old_prefix = 'read_count (' + str(labels_df["PrepLabel"][i])
        old_name = full_df.columns[full_df.columns.str.startswith(old_prefix)][0] #Find the old name from the prefix
        if labels_df["ID"][i][0:2] == 'In':
            tuples.append(('0d','In'))
        else:
            tuples.append(('6d',labels_df["ID"][i][0:-2]))
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
