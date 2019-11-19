from file_parser import file_parser
from gene_parser import gene_parser
from graph_filters import graph_filters
from resnik_calculator import resnik_calculator

import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from multiprocessing import Pool, freeze_support
import os

def delete_node_by_path(node_list, graph, graph_nodes):
        to_remove = []
        checked = []
        counter = 0
        for t in node_list:
            for s in graph_nodes:
                if s in checked:
                    continue
                if(nx.has_path(graph, s, t) == False):
                    to_remove.append(s)
                else:
                    checked.append(s)
            counter += 1
            if(counter % 300 == 0):
                print(str(counter) + "/" + str(len(node_list)))
        to_remove = list(set(graph_nodes) - set(checked))
        #print("Total irrelevant nodes removed is: %u" %len(to_remove))
        return to_remove    

G2 = nx.DiGraph()
G = nx.DiGraph()
fp = file_parser("./Data/Master_Go.obo")
gp = gene_parser("./Data/sgd_no_header.txt")

all_go_nodes = fp.generate_nodes_names()
gene_go_nodes = gp.create_go_nodes_from_gene()
valid_go_nodes = [n for n in gene_go_nodes if n in all_go_nodes]
print("Total valid nodes: %u" %len(valid_go_nodes))
valid_edges = gp.create_edge_list()
print("Total valid node-gene edge: %u" %len(valid_edges))
G2.add_edges_from(valid_edges)

fp = file_parser("./Data/Master_Go.obo")
nodes_names_list = fp.generate_nodes_names()
is_a, part_of, has_part, regulates = fp.generate_edges()
G.add_nodes_from(nodes_names_list, attribute = "GO_Term")
# If is_GO is 1, it means the edge is connecting two GO nodes. 
# If 0, it means the edge is connecting a Gene Node with an GO node
G.add_edges_from(is_a, is_GO = 1, is_Gene = 0, rel = "NEdge")
G.add_edges_from(part_of, is_GO = 1, is_Gene = 0, rel = "NEdge")
G.add_edges_from(has_part, is_GO = 1, is_Gene = 0, rel = "NEdge")
G.add_edges_from(regulates, is_GO = 1, is_Gene = 0, rel = "NEdge")

if __name__ == "__main__":
    print("in main")
    freeze_support()

    # MT for better performance
    thread_num = os.cpu_count()-1
    np_list = np.array(all_go_nodes)
    list_of_chucks = np.array_split(np_list, thread_num)
    iter_arg = [(valid_go_nodes, G, x.tolist()) for x in list_of_chucks]
    p = Pool()
    result = p.starmap(delete_node_by_path, iter_arg)
    flat_list = [item for sublist in result for item in sublist]
    flat_list = list(set(flat_list))
    print("total removed is %u" %len(flat_list))
    problems = [x for x in valid_go_nodes if x not in list(G.nodes)]
    p.close()
    p.join()
    G.remove_nodes_from(flat_list)

gene_edges = gp.create_edge_list()
gene_nodes = gp.create_gene_nodes_list()

G.add_nodes_from(gene_nodes, attribute = "Gene_Term")
G.add_edges_from(gene_edges, is_GO = 0, rel = "GDirect", is_Gene = 1)
gp.gene_propagation(G)

print("Graph Completed! \nTotal nodes: {} \nTotal edges: {}".format(len(G.nodes), len(G.edges)))

print("Applying filters ...")
gf = graph_filters(G, gp)
gf.build_GO_Gene_dict()
repeated_nodes_removed = gf.filter_repeated_nodes_2()
gene_few_nodes = gf.remove_nodes_with_too_few_genes_2(5)
gene_many_nodes = gf.remove_nodes_with_too_many_genes_2(0.3)

print("Filters Applied")

print("Generating file used for Resnik distance calculation...")
# generate GO nodes still in the graph
go_nodes = [x for x,d in G.nodes(data = True) if d['attribute'] == 'GO_Term']
with open ('termStats.txt', 'a') as f:
    for curr_term in go_nodes:
        curr_term_genes = nx.descendants(G, curr_term)
        #filter out go_terms
        curr_term_genes = [g for g in curr_term_genes if g in gene_nodes]
        total_gene_num = len(curr_term_genes)
        f.write("%s " %curr_term)
        f.write("%s " %total_gene_num)
        for g in curr_term_genes:
            f.write("%s," %g)
        f.write('\n')

print("File Generated")
