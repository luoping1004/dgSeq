import networkx as nx
G = nx.Graph()
dgDic = {}
fmap = open("D:/Data/OMIM/OMIM_v0.9.txt")
fmap.readline()

for line in fmap:
    data = line[:-1].split('\t')
    gene = data[1].split(', ')
    if data[0] in dgDic:
        dgDic[data[0]] = list(set(dgDic[data[0]]).union(set(gene)))
    else:        
        dgDic[data[0]] = gene

fmap.close()

Disease = list(dgDic.keys())
for disease in Disease:
    gene = dgDic[disease]
    for g in gene:
        G.add_edge(g,disease)

for i in range(len(Disease)):
    gene = dgDic[Disease[i]]
    for j in range(i+1,len(Disease)):
        for g in gene:
            if g in dgDic[Disease[j]]:
                G.add_edge(Disease[i],Disease[j])
                break


fdg = open("breast_cancer_gene_census.tsv")
fdg.readline()
dg = []
for line in fdg:
    data = line[:-1].split('\t')
    dg.append(data[0])

fdg.close()

GeneList = G.nodes()
for gene in dg:
    if gene in GeneList:
        GeneList.remove(gene)
for node in Disease:
    GeneList.remove(node)

collect = {}
nopath = []
for gene in GeneList:
    if nx.has_path(G,source=gene,target='114480'):#'114480' represents breast cancer
        length = nx.shortest_path_length(G,source=gene,target='114480')
        if length not in collect.keys():
            collect[length] = []
            collect[length].append(gene)
        else:
            collect[length].append(gene)
    else:
        nopath.append(gene)
print(collect.keys())

fdgn = open('dgn_shortest_path_norm.txt','w')
for k in collect.keys():
    gene = collect[k]
    fdgn.write('{}'.format(k))
    fdgn.write('\t')
    for i in range(len(gene)-1):
        fdgn.write(gene[i])
        fdgn.write(',')
    fdgn.write(gene[-1])
    fdgn.write('\n')

fdgn.close()

fnop = open('dgn_no_path_norm.txt','w')
for i in range(len(nopath)-1):
    fnop.write(nopath[i])
    fnop.write(',')
fnop.write(nopath[-1])
fnop.write('\n')
fnop.close()
















