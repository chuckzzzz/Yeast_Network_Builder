import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class gene_parser:
    def __init__(self, GO_Annotation_url):
        GO_data = pd.read_csv(GO_Annotation_url, sep = "\t", error_bad_lines=False)
        column_names = ["Database","Object ID","Object Symbol","Qualifier","GO ID","DB:Reference",
                "Evidence Code","With [or] From","Aspect","DB Object Name","DB Object Synonym",
                "DB Object Type","Taxon","Date","Assigned By","Annotation Extension","Gene Product Form ID"]
        GO_data.columns = column_names

        gene_list = GO_data["Object Symbol"]
        GO_ID_list = GO_data["GO ID"]
        self.genes = gene_list
        self.go_dict = {}
        
        for gene, ID in zip(gene_list, GO_ID_list):
            if ID not in self.go_dict:
                self.go_dict[ID] = []
                self.go_dict[ID].append(gene)
            else:
                cur_list = self.go_dict[ID]
                if gene not in cur_list:
                    cur_list.append(gene)
        # add first element 
        if("SPT23" not in self.go_dict["GO:0003674"]):
            self.go_dict["GO:0003674"].append("SPT23")
    
        
    def create_edge_list(self):
        result = [] #node, gene
        for key, value in self.go_dict.items():
            curr_GO_id = key
            curr_gene_list = value
            for gene in curr_gene_list:
                result.append((curr_GO_id, gene))
        result = list(set(result))
        return result
    
    def create_go_nodes_from_gene(self):
        nodes = self.go_dict.keys()
        return nodes
    
    def create_gene_nodes_list(self):
        gene_list = self.genes.tolist()
        unique_genes = list(set(gene_list))
        return unique_genes
    
    def gene_propagation(self, G):
        counter = 0
        unique_genes = self.create_gene_nodes_list()
        for gene in unique_genes:
            counter += 1
            if(counter %300 ==0):
                print(str(counter)+"/6300")
            curr_gene_anc = nx.ancestors(G, gene)
            curr_gene_par = G.predecessors(gene)
            indirect_par = curr_gene_anc.difference(curr_gene_par)
            #remove the direct parents
            #add the rest as indirect parents
            for a in indirect_par:    
                G.add_edge(a, gene, is_GO = 0, rel = "GIndirect", is_Gene = 1)
    