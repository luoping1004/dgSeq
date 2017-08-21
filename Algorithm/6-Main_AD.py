# This program compute the AUC of dgSeq on Alzheimer's Disease (AD) and output the
# average AUC value and a numpy matrix which contains the probabilities of
# all cross validation genes being disease-associated.

import numpy as np
import copy
from func_exclude_gene import ExcludeGene
from func_feature_generate import GetFeature
from func_minimum import FindMin
from func_label_graph import LabelGraph
from sklearn.metrics import roc_auc_score

#fdg = open("Alzheimer_dg_union.txt")
#fdg = open("50+omim_union.txt")
#fdg = open("182_Alz_genes.txt")
fdg = open("dg_top50.txt")
dg = fdg.readline()[:-1].split(',')
fdg.close()



fnode = open("GeneList_Alz.txt")
GeneList = []
for line in fnode:
    GeneList.append(line[:-1])
fnode.close()

dg = ExcludeGene(dg,GeneList)

#print(len(dg))

#none disease gene based on shortest path in normal disease-gene network

fdgnsp = open("dgn_shortest_path_norm_50.txt")
#dgn1norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
dgn2norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
dgn3norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
#dgn4norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
#dgn5norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
#dgn6norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
fdgnsp.close()

#fdgnnp = open("C:\\Users\\pil562\\Downloads\\Python\\BIBM2016_special_issue\\TCGA-THCA\\dgn_no_path_norm.txt")
fdgnnp = open("dgn_no_path_norm_50.txt")
dgnnpnorm = fdgnnp.readline()[:-1].split(',')
fdgnnp.close()

#none disease gene based on shortest path in differential network

fdgndiff = open("dgn_shortest_path_diff_50.txt")
dgn3diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
dgn4diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
dgn5diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
#dgn3diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
dgn100diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
fdgndiff.close()

dgnnorm = dgnnpnorm

dgnnorm = ExcludeGene(dgnnorm,GeneList)
#save the dgn
with open("AD_dgn.txt",'w') as fnw:
    for i in range(len(dgnnorm)-1):
        fnw.write(dgnnorm[i])
        fnw.write(',')
    fnw.write(dgnnorm[-1])
    fnw.write('\n')
    

dgndiff = dgn4diff+dgn5diff+dgn100diff


dgnpotential = list(set(dgnnorm).intersection(dgndiff))
alldgn = list(set(dgnnorm).union(set(dgndiff)))

print(len(dgnpotential))

wMatrix = np.load("Weight_Matrix_Pearson.npy")
Nn = 5
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

import networkx as nx
Gnorm = nx.Graph()
fppi = open("InBio_Map_Alz.txt")
for line in fppi:
    data = line[:-1].split('\t')
    Gnorm.add_edge(data[0],data[1])
fppi.close()

Gdiff = nx.Graph()
Gdiff.add_nodes_from(GeneList)
for i in range(m):
    for j in range(i+1,m):
        if newMatrix[i,j] != 0:
            Gdiff.add_edge(GeneList[i],GeneList[j],weight=newMatrix[i,j])
            
#print(Gnorm.number_of_nodes())
#print(Gdiff.number_of_nodes())
#print(len(dg))
#print(len(alldgn))


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
np.save("prob_AD.npy",probS)