#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd


in_fold = 'output/'
out_fold = 'output/'
padj_thresh = 0.05
lfc_thresh = 0.5
run_plots = True

#INH
fname = in_fold+'intersection_INH.csv'

df = pd.read_csv(fname)[['uid','H+++_7d','H+++_14d']]
df = df.rename({'H+++_7d':'LFC1','H+++_14d':'LFC2'},axis='columns')

plt.subplots()
df = df[(df['LFC1'].abs() > lfc_thresh) & (df['LFC2'].abs() > lfc_thresh)]
df = df.sort_values(by='LFC2',ascending=True)
print('Number of significant hypersusceptible mutants (INH):',sum(df['LFC2'] < 0))
print('Number of significant hypertolerant mutants (INH):',sum(df['LFC2'] > 0))
plt.bar(df['uid'],df['LFC2'],1,edgecolor='k',color='g')
plt.xlim(-0.5,5)
plt.ylim(-6,6)
plt.xlabel('Disrupted Genes')
plt.ylabel('Log$_2$(Fold-Change) (day 14, Highest Dose)')
plt.axhline(y=0,color='k')
plt.savefig(out_fold+'INH_PBS_barplot.png',bbox_inches='tight')

#RMP
fname = in_fold+'intersection_RMP.csv'

df = pd.read_csv(fname)[['uid','R+++_7d','R+++_14d']]
df = df.rename({'R+++_7d':'LFC1','R+++_14d':'LFC2'},axis='columns')

plt.subplots()
df = df[(df['LFC1'].abs() > lfc_thresh) & (df['LFC2'].abs() > lfc_thresh)]
df = df.sort_values(by='LFC2',ascending=True)
print('Number of significant hypersusceptible mutants (RMP):',sum(df['LFC2'] < 0))
print('Number of significant hypertolerant mutants (RMP):',sum(df['LFC2'] > 0))
print(df.iloc[0:10])
#print(df.iloc[-11:-1])
plt.bar(df['uid'],df['LFC2'],1,edgecolor=None,color='r')
#plt.xlim(-0.5,5)
plt.ylim(-6,6)
plt.xlabel('Disrupted Genes')
plt.ylabel('Log$_2$(Fold-Change) (day 14, Highest Dose)')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.axhline(y=0,color='k')
plt.savefig(out_fold+'RMP_PBS_barplot.png',bbox_inches='tight')
plt.show()

