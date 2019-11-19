import networkx as nx 
import matplotlib 
import numpy
class file_parser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_list = None
    
    def file_to_list(self):
        self.file_list = [line.rstrip('\n') for line in open(self.file_path)]
    
    def parse_GO_ID(self, line):
        result = line.split()
        try:
            return result[1]
        except IndexError:
            raise Exception("Error: " + " ".join(line))
    
    def generate_nodes_names(self):
        self.file_to_list()
        result = []
        sep = "[Term]"
        sep_met = False
        for line in self.file_list:
            if(sep in line):
                sep_met = True
                continue
            if(sep_met == True):
                node_name = self.parse_GO_ID(line)
                result.append(node_name)
                sep_met = False
        print("Total nodes: %u " %len(result))
        return result


G = nx.DiGraph()
fp = file_parser("./Data/Master_Go.obo")
nodes_names_list = fp.generate_nodes_names()
G.add_nodes_from(nodes_names_list, attribute = "GO_Term")
nx.draw(G, with_labels = True)