#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd


in_fold = 'output/'
out_fold = 'output/'
padj_thresh = 0.05
lfc_thresh = 0.5
run_plots = True

#INH
fname = in_fold+'summary_6d_INH.csv'

df = pd.read_csv(fname)[['uid','H+','pval','pval-adj (BH)']]
df = df.rename({'pval-adj (BH)':'padj','H+':'LFC'},axis='columns')

plt.subplots()
df = df[((df['padj'] < padj_thresh) & (df['LFC'].abs() > lfc_thresh))]
df = df.sort_values(by='LFC',ascending=True)
print('Number of significant hypersusceptible mutants (INH):',sum(df['LFC'] < 0))
print('Number of significant hypertolerant mutants (INH):',sum(df['LFC'] > 0))
plt.bar(df['uid'],df['LFC'],1,edgecolor='k',color='g')
plt.xlim(-0.5,5)
plt.ylim(-6,6)
plt.xlabel('Disrupted Genes')
plt.ylabel('Log$_2$(Fold-Change)')
plt.axhline(y=0,color='k')
plt.savefig(out_fold+'INH_barplot.png',bbox_inches='tight')

#RMP
fname = in_fold+'summary_6d_RMP.csv'

df = pd.read_csv(fname)[['uid','R+','pval','pval-adj (BH)']]
df = df.rename({'pval-adj (BH)':'padj','R+':'LFC'},axis='columns')

plt.subplots()
df = df[((df['padj'] < padj_thresh) & (df['LFC'].abs() > lfc_thresh))]
df = df.sort_values(by='LFC',ascending=True)
print('Number of significant hypersusceptible mutants (RMP):',sum(df['LFC'] < 0))
print('Number of significant hypertolerant mutants (RMP):',sum(df['LFC'] > 0))
print(df.iloc[0:10])
#print(df.iloc[-11:-1])
plt.bar(df['uid'],df['LFC'],1,edgecolor=None,color='r')
#plt.xlim(-0.5,5)
plt.ylim(-6,6)
plt.xlabel('Disrupted Genes')
plt.ylabel('Log$_2$(Fold-Change)')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.axhline(y=0,color='k')
plt.savefig(out_fold+'RMP_barplot.png',bbox_inches='tight')
plt.show()

