import matplotlib.pyplot as plt
import pandas as pd

import matplotlib
import seaborn as sns
import scipy

import sys
sys.setrecursionlimit(100000)

CCC_BNEW = pd.read_csv('/home/users/testme/My_/DIGIST/cell_type_CCC_vs_BNEW.csv')
CRC_BNEW = pd.read_csv('/home/users/testme/My_/DIGIST/cell_type_CRC_vs_BNEW.csv')
CRT_BNEW = pd.read_csv('/home/users/testme/My_/DIGIST/cell_type_CRT_vs_BNEW.csv')

gene_union = set()
gene_intersection = set()

for df in [CCC_BNEW, CRC_BNEW, CRT_BNEW]:
    df_ = df.loc[~df['padj'].isna()]
    df_ = df_[df_['padj'] < 0.05]
    gene = df_['Unnamed: 0']
    gene_union = gene_union.union(gene)
    if len(gene_intersection) == 0:
        gene_intersection = set(gene)
    gene_intersection = gene_intersection.intersection(gene)
    
cts = pd.read_csv('/home/users/testme/My_/DIGIST/merged_count_matrix_normalized.csv')
cts.index = cts['Unnamed: 0']
cts.index.name = 'gene'

cts_union = cts[cts['Unnamed: 0'].isin(gene_union)].drop(columns=['Unnamed: 0'], axis=1)
cts_inter = cts[cts['Unnamed: 0'].isin(gene_intersection)].drop(columns=['Unnamed: 0'], axis=1) # it is not used, but can replace cts_union with cts_inter in the following code

g = sns.clustermap(cts_union, 
               cmap = matplotlib.cm.get_cmap('Spectral_r'), 
               row_cluster=True, 
               col_cluster=False, 
               metric = 'correlation',
               )
plt.title('union')

den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_row.linkage)
scipy.cluster.hierarchy.fcluster(g.dendrogram_row.linkage, t=1.2, criterion='distance')

g = sns.clustermap(cts_union_clust.drop(columns=['clust']), 
               cmap = matplotlib.cm.get_cmap('Spectral_r'), 
               row_cluster=False, 
               col_cluster=False, 
               row_colors=cts_union_clust['clust'].map({1 : 'r', 2 : 'g'})
               )
plt.title('union')

cts_ = cts_union_clust[cts_union_clust['clust'] == 1]
g_ = sns.clustermap(cts_.drop(columns=['clust'], axis=1), 
            cmap = matplotlib.cm.get_cmap('Spectral_r'), 
            row_cluster=True, 
            col_cluster=False,
)
plt.title(f'union_{1}')

den = scipy.cluster.hierarchy.dendrogram(g_.dendrogram_row.linkage)
array = scipy.cluster.hierarchy.fcluster(g_.dendrogram_row.linkage, t=33, criterion='distance') # The cutoff t=33 follows from the interative search pd.Series(array).value_counts()

cts_['clust'] = array 
cts_union_clust_1 = cts_.sort_values(by = 'clust')

temp = cts_union_clust_1['clust'].map({11 : 'r', 6 : 'g', 7 : 'b'})
temp[temp.isna()] = 'c'

g = sns.clustermap(cts_union_clust_1.drop(columns=['clust']), 
               cmap = matplotlib.cm.get_cmap('Spectral_r'), 
               row_cluster=False, 
               col_cluster=False, 
               row_colors=temp
               )
plt.title('union_1_')

cts_ = cts_union_clust[cts_union_clust['clust'] == 2]
g_ = sns.clustermap(cts_.drop(columns=['clust'], axis=1), 
            cmap = matplotlib.cm.get_cmap('Spectral_r'), 
            row_cluster=True, 
            col_cluster=False,
)
plt.title(f'union_{2}')

den = scipy.cluster.hierarchy.dendrogram(g_.dendrogram_row.linkage)
array = scipy.cluster.hierarchy.fcluster(g_.dendrogram_row.linkage, t=33, criterion='distance')

cts_['clust'] = array 
cts_union_clust_2 = cts_.sort_values(by = 'clust')

temp = cts_union_clust_2['clust'].map({9 : 'r', 20 : 'g', 18 : 'b'})
temp[temp.isna()] = 'c'

g = sns.clustermap(cts_union_clust_2.drop(columns=['clust']), 
               cmap = matplotlib.cm.get_cmap('Spectral_r'), 
               row_cluster=False, 
               col_cluster=False, 
               row_colors=temp
               )
plt.title('union_2_')
