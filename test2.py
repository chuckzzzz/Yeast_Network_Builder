import numpy as np
import networkx as nx 

#        1
#    2       3

G = nx.DiGraph()
nodes = [1,2,3,4]
G.add_nodes_from(nodes)
edges = [(1,2),(1,3), (2,4)]
G.add_edges_from(edges)
print(G.nodes)

nodes = [1, 2]
with open ('test.txt', 'a') as f:
    for curr_term in nodes:
        curr_term_genes = nx.descendants(G, curr_term)
        total_gene_num = len(curr_term_genes)
        f.write("%s " %curr_term)
        f.write("%s " %total_gene_num)
        for g in curr_term_genes:
            f.write("%s," %g)
        f.write('\n')
