import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gene_parser import gene_parser


class graph_filters:
    def __init__(self, graph, gene_parser):
        self.graph = graph
        self.go_nodes_names = [x for x,y in self.graph.nodes(data = True) if y['attribute'] == "GO_Term"]
        self.unique_genes = gene_parser.create_gene_nodes_list()
        self.go_dict = {}
    
    def get_leaf_nodes_list(self):
        result = [x for x in self.graph.nodes if self.graph.out_degree(x)==0]
        return result

    # all nodes with less than threshold number (int) of genes will be abandoned
    def filter_nodes_by_genes(self, threshold):
        counter = 0
        to_remove = []
        for node_id in self.go_nodes_names:
            # Get all recursive children of the current node
            cur_des = nx.descendants(self.graph, node_id)
            curr_gene_count = 0
            # Count the number of children that's a gene
            for d in list(cur_des):
                cur_node = self.graph.node[d]
                if (cur_node["attribute"] == "Gene_Term"):
                    curr_gene_count += 1
            # Check if the number of gene children is smaller than the threshold parameter
            if(curr_gene_count < threshold):
                # attach deleted node's gene to ancestors
                cur_pre = self.graph.predecessors(node_id)
                cur_sec = self.graph.successors(node_id)
                # add all children of the current node to all parents of the current node 
                for p in cur_pre:
                    for s in cur_sec:
                        self.graph.add_edge(p, s, attribute = "GO_Gene", is_GO = 0, rel = " ")
                self.graph.remove_node(node_id)
                to_remove.append(node_id)
                counter += 1
        for r in to_remove:
            if r in self.go_nodes_names:    
                self.go_nodes_names.remove(r)
        print("Total remove by too few genes is %u" %counter)
    
    def subset_nodes_by_attr(self, attr):
        nodes = [x for x,y in self.graph.nodes(data=True) if y['attribute'] == attr]
        return nodes
# while true loop, inside for loop 
# father only indirect edge, only one child, child only one parent
# add child to all parent ancestors 
# judge single parent single children(G) return true 
    def filter_repeated_nodes(self):
        counter = 0
        for node_id in self.go_nodes_names:
            cur_GO_out = self.graph.out_degree(node_id, weight = "is_GO")
            # there is only one outdegree for GO_Node
            if(cur_GO_out == 1):
                # loop through its successors
                cur_suc = self.graph.successors(node_id)
                for suc in list(cur_suc).copy():
                    s_node = self.graph.node[suc]
                    # if the current successor is a GO_term, which there should only be one, check if its in_degree is 1
                    if(s_node["attribute"] == "GO_Term"):
                        cur_suc_in_degree = self.graph.in_degree(suc, weight = "is_GO")
                        # if the in degree is also 1, add all of its gene edges to the original node and delete this node
                        if(cur_suc_in_degree == 1):
                            counter += 1
                            # found repeat node, delete this successor and attach all of its successors to the node
                            suc_of_suc = list(self.graph.successors(suc))
                            new_gene_edges = [(node_id,y) for y in suc_of_suc if self.graph.nodes[y]["attribute"] == "Gene_Term"]
                            new_go_edges = [(node_id,y) for y in suc_of_suc if self.graph.nodes[y]["attribute"] == "GO_Term"]
                            self.graph.add_edges_from(new_gene_edges, is_GO = 0)
                            self.graph.add_edges_from(new_go_edges, is_GO = 1)
                            self.graph.remove_node(suc)
                            self.go_nodes_names.remove(suc)

        print("Total remove by repeated is %u" %counter)
        
        # genes with over threshold percent will be removed 
    def remove_nodes_with_too_many_genes(self, threshold):
        threshold_num = (int)(threshold * len(self.unique_genes))
        counter = 0
        print("Threshold: %u" %threshold_num)
        for node_id in self.go_nodes_names:
            cur_des = nx.descendants(self.graph, node_id)
            curr_gene_num = 0
            for des in list(cur_des):
                if(self.graph.nodes[des]["attribute"] == "Gene_Term"):
                    curr_gene_num += 1
            if(curr_gene_num > threshold_num):
                #cur_suc = list(self.graph.successors(node_id)
                #new_gene_edges = [(node_id,y) for y in cur_suc if self.graph.nodes[y]["attribute"] == "Gene_Term"]
                self.graph.remove_node(node_id)
                self.go_nodes_names.remove(node_id)
                counter += 1
        print("Total remove by too many genes is %u" %counter)
    
    #use this method with G'
    def build_GO_Gene_dict(self):
        for node_id in self.graph.nodes:
            if(self.graph.nodes[node_id]["attribute"] == "Gene_Term"):
                continue
            cur_gene_count = self.graph.out_degree(node_id, weight = "is_Gene")
            self.go_dict[node_id] = cur_gene_count
        return None
    
     # while (Found == True)
            # loop through the nodes
            # check node out_degree
            # when out_degree == 1
                # check child in_degree
                # when child in_degree ==1
                    # Found == True
                    # get parent node parents
                    # create edge to all grand parents node
                    # remove the parent node 
                    # break
            # Found = False 
    def filter_repeated_nodes_2(self):
        found = True
        first_loop = True
        counter = 0
        node_id_to_remove = None
        result = []
        while(found):
            if(first_loop == False):
                self.graph.remove_node(node_id_to_remove)
                del self.go_dict[node_id_to_remove]
                result.append(node_id_to_remove)
            first_loop = False
            for node_id in self.graph.nodes:
                cur_out_degree = self.graph.out_degree(node_id, weight = "is_GO")
                if(cur_out_degree == 1):
                    cur_suc = self.graph.successors(node_id)
                    for s in cur_suc:
                        if self.graph.nodes[s]['attribute'] == "GO_Term":
                            cur_child_in_degree = self.graph.in_degree(s, weight = "is_GO")
                            if (cur_child_in_degree == 1):
                                cur_grand_parents = self.graph.predecessors(node_id)
                                for gp in cur_grand_parents:
                                    self.graph.add_edge(gp, s, is_GO = 1, is_Gene = 0, rel = "NEdge" )
                                node_id_to_remove = node_id
                                counter += 1
                                break
            found = False
        print("Total repeated nodes removed is %u" %counter)
        return result
    
    def remove_nodes_with_too_many_genes_2(self, threshold):
        threshold = int(threshold * len(self.unique_genes))
        to_remove = []
        for node_id in self.graph.nodes:
            if node_id not in self.go_dict.keys():
                continue
            curr_node_gene_num = self.go_dict[node_id]
            if(curr_node_gene_num > threshold):
                to_remove.append(node_id)                
                del self.go_dict[node_id]
        print("Updating edges bofore removal...")
        self.update_edges_before_removal(to_remove, "rel")
        self.graph.remove_nodes_from(to_remove)    
        # update edge attribute with all the parents 
        print("Total nodes removed by too many genes is %u" %len(to_remove))
        
        return to_remove
    
    def check_for_repeated(self):
        for node_id in self.graph.nodes:
            cur_out_degree = self.graph.out_degree(node_id, weight = "is_GO")
            if(cur_out_degree == 1):
                    cur_suc = self.graph.successors(node_id)
                    for s in cur_suc:
                        if self.graph.nodes[s]['attribute'] == "GO_Term":
                            cur_child_in_degree = self.graph.in_degree(s, weight = "is_GO")
                            if (cur_child_in_degree == 1):
                                print(node_id, s)
                                return True
        return False
            
            
            
    def remove_nodes_with_too_few_genes_2(self, threshold):
        to_remove = []
        for node_id in self.graph.nodes:
            if node_id not in self.go_dict.keys():
                continue
            curr_node_gene_num = self.go_dict[node_id]
            if(curr_node_gene_num < threshold):
                to_remove.append(node_id)
                del self.go_dict[node_id]
        print("Updating edges bofore removal...")
        self.update_edges_before_removal(to_remove, "rel")
        self.graph.remove_nodes_from(to_remove)    
        # update edge attribute with all the parents 
        print("Total nodes removed by too many genes is %u" %len(to_remove))
        
        return to_remove
    
    def get_children_by_attr(self, node_id, attr, edge_attr = None):
        suc = self.graph.successors(node_id)
        result = [x for x in suc if self.graph.nodes[x]['attribute'] == attr]
        if(edge_attr == None):
            return result
        else:
            edges = [(node_id, s) for s in result]
            result_e = set([s for (node_id, s) in edges if self.graph.edges[(node_id, s)]['rel'] == edge_attr]) 
            return result_e
    
    def update_edges_before_removal(self, node_list, edge_attr): #edge_attr = "rel"
        total = len(node_list)
        print(str(total))
        counter = 0
        for node_id in node_list:
            gene_list = self.get_children_by_attr(node_id, "GO_Term", "GDirect")
            cur_parent_list = self.graph.predecessors(node_id)
            if(len(gene_list)> 0):
                print(node_id + ":" + str(len(gene_list)))
            for g in gene_list:
                for p in cur_parent_list:
                    #print("current parent: " + p)
                    p_direct_children = self.graph.successors(p)
                    found = False
                    for c in p_direct_children:
                        #print("current direct child: " + c)
                        if c == node_id:
                            continue
                        if(g in nx.descendants(self.graph, c)):
                            found = True
                            break
                    if(found == False):
                        self.graph[p][g]["rel"] = "GDirect"
            counter += 1
            if(counter % 1000 == 0):
                print(str(counter) + "/" + str(total))
        return None