# This program compute the AUC of dgSeq on Breast Cancer (BC) and output the
# average AUC value and a numpy matrix which contains the probabilities of
# all cross validation genes being disease-associated.


import numpy as np
import networkx as nx
import copy
from func_exclude_gene import ExcludeGene
from func_feature_generate import GetFeature
from func_minimum import FindMin
from func_label_graph import LabelGraph
from sklearn.metrics import roc_auc_score

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
with open("dgn_shortest_path_norm_bc.txt") as fdgnsp:
    dgn1norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn2norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn3norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn4norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn5norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn6norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
    dgn7norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')

# This file contains none disease genes in disease-gene network which are not
# connected to the target disease
with open("dgn_no_path_norm_bc.txt") as fdgnnp:
    dgnnpnorm = fdgnnp.readline()[:-1].split(',')

# This file contains none disease genes based on the differential network
# The number in each line also denotes the length of the shortest path
with open("dgn_shortest_path_diff_bc.txt") as fdgndiff:
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

# The codes in the following comment section import the PPI network and the generate the differential network. Since the weight matrix P
# is too large to upload, we directly load the resulting PPI networks into the program.
"""
# The correlation difference matrix P
wMatrix = np.load("Weight_Matrix_BC.npy")
Nn = 5 # Top 5 neighbors are chosen to build the differential network
m,m = wMatrix.shape
newMatrix = np.zeros((m,m))
for i in range(m):
    test = wMatrix[i]
    temp = np.argpartition(-test, Nn)
    result_args = temp[:Nn]
    temp = np.partition(-test, Nn)
    result = -temp[:Nn]
    for j in range(Nn):
        newMatrix[i,result_args[j]] = result[j]
        newMatrix[result_args[j],i] = result[j]


# The PPI network
Gnorm = nx.Graph()
with open("InBio_Map_ppi_BC.txt") as fppi:
    for line in fppi:
        data = line[:-1].split('\t')
        Gnorm.add_edge(data[0],data[1])
# The differential network
Gdiff = nx.Graph()
for i in range(m):
    for j in range(i+1,m):
        if newMatrix[i,j] != 0:
            Gdiff.add_edge(GeneList[i],GeneList[j],weight=newMatrix[i,j])
"""

Gnorm = nx.read_gml("ppi_bc.gml.gz")
Gdiff = nx.read_gml("dif_bc.gml.gz")

# Compute the average AUC
Round = 100
auc100 = np.zeros(Round)
dgLength = len(dg)
prob = np.zeros(2*dgLength)
probS = np.zeros((Round,2*dgLength))
dg_tmp = copy.deepcopy(dg)
y = np.concatenate((np.ones(dgLength),np.zeros(dgLength)))
AllLength = 2*dgLength

import time
start = time.clock()

for r in range(Round):
    print(r)    
    # Randomly select len(dg) non-disease genes from benchmark set
    dgn = []
    while len(dgn) < dgLength:
        temp = np.random.choice(dgnpotential)
        if temp not in dgn:
            dgn.append(temp)
    AllGene = dg+dgn
       
    # leave one out
    for i in range(AllLength):      
        LOO = AllGene[i]
        
        # Label the graph
        if i < dgLength:
            dg_tmp.remove(LOO)
            LabelGraph(Gnorm, Gdiff, dg_tmp, alldgn)
            dg_tmp.append(LOO)
        else:
            alldgn.remove(LOO)
            LabelGraph(Gnorm, Gdiff, dg_tmp, alldgn)
            alldgn.append(LOO)
            
        # Extract feature
        X = GetFeature(Gnorm, Gdiff, AllGene, AllLength)
        indices = np.arange(2*dgLength)
        train_index = indices[np.logical_not(indices==i)]
        X_train, y_train = X[train_index], y[train_index]
        X_test = X[i]
        
        # Compute parameters
        res = FindMin(X_train,y_train)
        FW = np.dot(X_test,res.x)
        prob[i] = np.exp(FW)/(1+np.exp(FW))
     
    auc100[r] = roc_auc_score(y, prob)
    probS[r] = prob

end = time.clock()
print ('run:',(end - start)/60,'min')

#save the probabilities for drawing the ROC curve
print(np.mean(auc100))
np.save("prob_BC.npy",probS)
