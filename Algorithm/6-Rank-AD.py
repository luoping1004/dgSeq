import numpy as np
import networkx as nx

from func_exclude_gene import ExcludeGene
from func_feature_generate import GetFeature
from func_minimum import FindMin
from func_label_graph import LabelGraph

fdg = open("Alzheimer's_disease_genes.txt")
dg = fdg.readline()[:-1].split(',')
fdg.close()



fnode = open("GeneList_Alz.txt")
GeneList = []
for line in fnode:
    GeneList.append(line[:-1])
fnode.close()

dg = ExcludeGene(dg,GeneList)

#none disease gene based on shortest path in normal disease-gene network

fdgnsp = open("dgn_shortest_path_norm_ad.txt")
dgn2norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
dgn3norm = (fdgnsp.readline()[:-1].split('\t')[1]).split(',')
fdgnsp.close()

fdgnnp = open("dgn_no_path_norm_ad.txt")
dgnnpnorm = fdgnnp.readline()[:-1].split(',')
fdgnnp.close()

#none disease gene based on shortest path in differential network
fdgndiff = open("dgn_shortest_path_diff_ad.txt")
dgn3diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
dgn4diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
dgn5diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
dgn100diff = (fdgndiff.readline()[:-1].split('\t')[1]).split(',')
fdgndiff.close()

dgnnorm = dgnnpnorm

dgnnorm = ExcludeGene(dgnnorm,GeneList)

dgndiff = dgn4diff+dgn5diff+dgn100diff

dgnpotential = list(set(dgnnorm).intersection(dgndiff))
alldgn = list(set(dgnnorm).union(set(dgndiff)))


Gnorm = nx.read_gml("ppi_ad.gml.gz")
Gdiff = nx.read_gml("dif_ad.gml.gz")

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

fuw = open("ad_unknown_gene_list.rnk",'w')
for i in range(100):
    fuw.write(unknownRe[i][0])
    fuw.write('\t')
    fuw.write('{}'.format(unknownRe[i][1]))
    fuw.write('\n')
fuw.close()

