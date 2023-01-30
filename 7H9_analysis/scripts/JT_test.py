#!/usr/bin/env python3
import time
import numpy as np
from itertools import combinations
from statsmodels.distributions.empirical_distribution import ECDF
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd

def main():
    #TODO: Test to ensure p-values are accurate
    #data_ls = [1,2,0,4,0,6,0,8,9,10,11,12]
    #group_ls = ['a','a','a','b','b','b','c','c','c','d','d','d']
    #expected_order = ['a','b','c','d']

    data_ls = [74,58,68,60,69,70,72,75,80,71,73,78,88,85,76]
    group_ls = ['a']*5 + ['b']*5 + ['c']*5
    expected_order = ['a','b','c']
    
    df = pd.DataFrame([data_ls]*100,columns = group_ls)

    t0 = time.time()
    p = JT_test_pd(df,expected_order,direction='lessthan',random_smooth=False)
    print('pvalue:',p)
    t1 = time.time()
    print(t1-t0)

def JT_test_pd(df,expected_order,direction='lessthan',random_smooth=False):
    #df is a pandas dataframe. Column names are used as the group names.
    #expected_order is a list containing the expected order of the group names (column names). 
    #direction is a string containing either "greaterthan" or "lessthan". If direction == "lessthan" then the alternative hypothesis expects that expected_order[0] < expected_order[1] < expected_order[2]...
    #Set random_smooth = True if you want to convert the discrete pvalue to a smooth uniformly distributed pvalue. This can be helpful when follow-on methods assume uniformly distributed p-values.
    #Output: Pandas series containing pvalues.
    if direction=='lessthan':
        pass
    elif direction=='greaterthan':
        expected_order = expected_order[::-1]
    else:
        exit('Direction: '+direction+' not recognized.')

    group_ls = list(df.columns)
    num_samp = len(group_ls)

    #Calculate ecdf for all possible ties at rank 1 through num_samp
    ecdfr_ls = [0]*(num_samp+1)
    ecdfl_ls = [0]*(num_samp+1)
    eval_plots = np.arange(0,50,0.5)
    for num_ties in range(num_samp+1):
        ecdfr_ls[num_ties],ecdfl_ls[num_ties] = JT_test_cdf([-1]*num_ties+list(range(num_samp-num_ties)),group_ls) #Note that first two cdfs will be identical - can't have "one" tied sample.

    #Calculate p-values using ecdf
    def get_pvalue(row,random_smooth=False):
        num_ties = len(row) - len(set(row))

        A = defaultdict(list)
        for k,v in zip(group_ls,row):
            A[k].append(v)
        Bprime = JT_statistic([A[k] for k in expected_order]) #Value of JT statistic for input data

        if random_smooth:
            pval = 1 - np.random.uniform(ecdfl_ls[num_ties](Bprime),ecdfr_ls[num_ties](Bprime))
        else:
            pval = 1 - ecdfl_ls[num_ties](Bprime)
        return pval

    pval_sr = df.apply(get_pvalue,axis=1,args=(random_smooth,))

    return pval_sr

def JT_test_cdf(data_ls,group_ls):
    #data_ls is a list containing values from samples
    #group_ls is a list containing group labels corresponding to each sample (ie if sample 1 is from group 'A' and was measured to be 0.7 than data_ls[0] = 0.7 and group_ls = 'A'. Values in group_ls must be allowable as keys in a dictionary.
    #Returns: ecdf as computed with statsmodels

    groups = list(set(group_ls))
    A = defaultdict(list)
    for k,v in zip(group_ls,data_ls):
        A[k].append(v)
    num_group = [0]*len(groups)
    for i,k in enumerate(groups):
        num_group[i] = len(A[k])

    #Generate All Random Combinations
    ind = list(range(len(data_ls))) #These are the indices of each element
    B = []
    data_ls_np = np.array(data_ls)
    for combo_ind in recursive_combination(ind,num_group):
        combo = [data_ls_np[k_ls] for k_ls in combo_ind]
        B.append(JT_statistic(combo))

    #Calculate ecdf
    #The JT distribution is fixed for all lists of unique numbers. If there are ties then the rank of the tied value(s) makes a difference. However for two lists of values with ties at the same ranks the distributions will be equal.
    ecdf_r = ECDF(B,side='right')
    ecdf_l = ECDF(B,side='left')
    return (ecdf_r,ecdf_l)

def recursive_combination(int_ls,num_ls):
    #Find all unique assignments of integers in int_ls to groups of size num_ls.
    #Assumes that all elements of int_ls are unique
    if len(int_ls) == num_ls[0]:
        yield [int_ls] #There is only one combination of values of int_ls
    else:
        for part_combo in combinations(int_ls,num_ls[0]):
            rem = list(set(int_ls) - set(part_combo))
            rem_lls_gen = recursive_combination(rem,num_ls[1:])
            for r in rem_lls_gen:
                lls = [list(part_combo)] + r
                yield lls

def JT_test_perm(data_ls,group_ls,N_perm,direction='lessthan'):
    #Input: data_ls is a list of values, group_ls is the group that each value belongs to, N_perm is the number of permutations that should be done.
    if direction=='lessthan':
        pass
    elif direction=='greaterthan':
        data_ls = data_ls[::-1]
        group_ls = group_ls[::-1]
    else:
        exit('Direction: '+direction+' not recognized.')

    #Check that group_ls is sorted
    #Calculate JT statistic (B)
    B = JT_statistic(data_ls,group_ls)
    #Permute samples, build list of JT statistic
    I_perm_ls = [None]*N_perm
    for i in range(N_perm):
        I_perm_ls[i] = (JT_statistic(np.random.permutation(data_ls),group_ls) >= B)

    return sum(I_perm_ls)/N_perm

def JT_statistic(lls):
    #Input: lls is a list of list. Each sublist is a group.
    #Note: This technically implements Terpstra's statistic - though there is a linear transform to get Jonckheere's statistics. This won't affect the results of the permutation test however.

    #Output: Calculated value of JT statistic
    #Calculate JT statistic
    k = len(lls)
    B = 0
    for i in range(k-1):
        for j in range(i+1,k):
            for r in range(len(lls[i])):
                for s in range(len(lls[j])):
                    if lls[i][r] < lls[j][s]:
                        D = 1
                    elif lls[i][r] == lls[j][s]:
                        D = 0.5
                    else:
                        D = 0
                    B = B + D
    return B

if __name__ == '__main__':
    main()
