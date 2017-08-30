import numpy as np
import networkx as nx

with open("breast_cancer_gene_census.tsv") as fdg:
    fdg.readline()
    dg = []
    for line in fdg:
        data = line.split('\t')
        dg.append(data[0])

with open("GeneList_BRCA.txt") as fnode:
    GeneList = []
    for line in fnode:
        data = line.split('\t')
        GeneList.append(data[1])

exclude = []
for gene in dg:
    if gene not in GeneList:
        exclude.append(gene)
for gene in exclude:
    dg.remove(gene)

Gdiff = nx.read_gml("dif_bc.gml.gz")
    
record = [None for col in range(2)]
record[1] = 100
collect = {}
for gene in GeneList:
    for dgene in dg:
        if nx.has_path(Gdiff,source=gene,target=dgene):
            length = nx.shortest_path_length(Gdiff,source=gene,target=dgene)
            if length < record[1]:
                record[0] = dgene
                record[1] = length
    if (record[1] > 2):
        if record[1] not in collect.keys():
            collect[record[1]] = []
            collect[record[1]].append(gene)
        else:
            collect[record[1]].append(gene)
    record[1] = 100
    
keyset = collect.keys()
keyset = sorted(keyset)


fdgn = open('dgn_shortest_path_diff_bc.txt','w')
for k in keyset:
    gene = collect[k]
    fdgn.write('{}'.format(k))
    fdgn.write('\t')
    for i in range(len(gene)-1):
        fdgn.write(gene[i])
        fdgn.write(',')
    fdgn.write(gene[-1])
    fdgn.write('\n')
fdgn.close()

