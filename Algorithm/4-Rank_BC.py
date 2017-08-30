import numpy as np
import networkx as nx

from func_exclude_gene import ExcludeGene
from func_feature_generate import GetFeature
from func_minimum import FindMin
from func_label_graph import LabelGraph

# Disease genes file downloaded from COSMIC
with open("breast_cancer_gene_census.tsv") as fdg:
    fdg.readline()
    dg = []
    for line in fdg:
        data = line.split('\t')
        dg.append(data[0])

# This file contins a list of genes generated through combining the original
# PPI network and BC RNA-Seq data. Only the protein coding genes with nonzero  
# expression values exist in this file.
with open("GeneList_BRCA.txt") as fnode:
    GeneList = []
    for line in fnode:
        data = line.split('\t')
        GeneList.append(data[1])


# This file contains none disease genes based on shortest path in disease-gene
# network (DGN). The number in each line denotes the length of the shortest path
with open("dgn_shortest_path_norm.txt") as fdgnsp:
    dgn1norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn2norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn3norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn4norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn5norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn6norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn7norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')

# This file contains none disease genes in disease-gene network which are not
# connected to the target disease
with open("dgn_no_path_norm.txt") as fdgnnp:
    dgnnpnorm = fdgnnp.readline()[:-1].split(',')

# This file contains none disease genes based on the differential network
# The number in each line also denotes the length of the shortest path
with open("dgn_shortest_path_diff.txt") as fdgndiff:
    dgn3diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
    dgn4diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
    dgn5diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
fdgndiff.close()

# dgnnorm stores non-disease genes collected from DGN
dgnnorm = dgn5norm+dgn6norm+dgn7norm+dgnnpnorm
# Genes not in the PPI network are removed
dgnnorm = ExcludeGene(dgnnorm,GeneList)

#dgndiff stores non-disease genes collected from the differential network
dgndiff = dgn4diff+dgn5diff

# Benchmark non-disease genes set
dgnpotential = list(set(dgnnorm).intersection(dgndiff))

# alldgn contains genes whose prior labels are zero
alldgn = list(set(dgnnorm).union(set(dgndiff)))

Gnorm = nx.read_gml("ppi_bc.gml.gz")
Gdiff = nx.read_gml("dif_bc.gml.gz")

# Compute the average AUC

import copy
unknown = copy.deepcopy(GeneList)

for gene in dg:
    unknown.remove(gene)
for gene in alldgn:
    unknown.remove(gene)

   
Round = 100
prob_unknown = np.zeros((len(unknown), Round))
dgLength = len(dg)
y = np.concatenate((np.ones(dgLength),np.zeros(dgLength)))

# Label the graph
LabelGraph(Gnorm, Gdiff, dg, alldgn)
# Extract feature for unknown genes
X_unknown = GetFeature(Gnorm, Gdiff, unknown, len(unknown))

import time
start = time.clock()

for r in range(Round):
    #print(r)    
    # Randomly select len(dg) non-disease genes from benchmark set
    dgn = []
    while len(dgn) < dgLength:
        temp = np.random.choice(dgnpotential)
        if temp not in dgn:
            dgn.append(temp)
    AllGene = dg+dgn
    AllLength = len(AllGene)
       
    # Extract feature
    X = GetFeature(Gnorm, Gdiff, AllGene, AllLength)
        
    # Compute parameters
    res = FindMin(X,y)
    
    # compute probabilities
    for i in range(len(unknown)):
        prob_unknown[i,r] = np.exp(np.dot(X_unknown[i],res.x))/(1+np.exp(np.dot(X_unknown[i],res.x)))

end = time.clock()
print ('run:',(end - start)/60,'min')

unknownDic = {}
for i in range(len(unknown)):
    unknownDic[unknown[i]] = np.mean(prob_unknown[i])
from operator import itemgetter
unknownRe = sorted(list(unknownDic.items()), key=itemgetter(1), reverse=True)

fuw = open("BC_top100.txt",'w')
for i in range(100):
    fuw.write(unknownRe[i][0])
    fuw.write('\n')
fuw.close()

