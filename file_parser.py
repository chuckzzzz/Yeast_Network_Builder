import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from multiprocessing import Pool
import os

class file_parser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.rels = ["is_a", "part_of", "has_part", "regulates"]
        self.black_list = ['def:', 'comment:']
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
                continue
        print("Total nodes: %u " %len(result))
        return result
    
    def parse_GO_id(self, data_id):
        result = data_id.split()
        if(len(result) == 0):
            return None
        try:
            return result[1]
        except IndexError:
            print("Error: "+data_id)
            raise Exception("Error: " + " ".join(data_id))
            
    def parse_relationship(self, rel):
        r_list = rel.split()
        if(len(r_list)>0):
            if r_list[0] in self.black_list:
                return (None, None)
        index = 0
        GO_id = None
        GO_rel = None
        for index in range (len(r_list)):
            if "GO:" in r_list[index]:
                GO_id = r_list[index]
            for rel in self.rels:
                if rel in r_list[index]:
                    GO_rel = rel
        return (GO_id, GO_rel)
    
    def delete_nodes_by_path(self, node_list, G):
        print(node_list)
        to_remove = []
        counter = 0
        for t in node_list:
            for s in G.nodes:
                if s in node_list:
                    continue
                if(nx.has_path(G, s, t) == False):
                    to_remove.append(s)
            counter += 1
            if(counter % 300 == 0):
                print(str(counter) + "/" + str(len(node_list)))
        to_remove = list(set(to_remove))                
        print("Total irrelevant nodes removed is: %u" %len(to_remove))
        return to_remove
    
    def delete_nodes_by_path_mt(self, node_list, G):
        start_time = time.time()
        thread_num = os.cpu_count()
        np_list = np.array(node_list)
        list_of_chucks = np.array_split(np_list, 6)
        iter_arg = [(x.tolist(), G) for x in list_of_chucks]
        p = Pool()
        result = p.map(self.delete_nodes_by_path, iter_arg)
        p.close()
        p.join()
        G.remove_nodes_from(result)
        end_time = str(time.time()-start_time)
        print("Total time used is {end_time}")
        

        
    def generate_edges(self):
        is_a = []
        part_of = []
        has_part = []
        regulates = []
        sep = "[Term]"
        stop = "[Typedef]"
        curr_id = None
        sep_met = False
        for curr in self.file_list:
            if curr == stop:
                break
            if sep_met == True:
                curr_id = self.parse_GO_id(curr)
                sep_met = False
                continue
            if curr == sep: 
                sep_met = True
                continue
            target_id, cur_rel = self.parse_relationship(curr)
            if cur_rel == None or target_id == None:
                continue
            if "is_a" in cur_rel:
                is_a.append((target_id, curr_id))
            if "part_of" in cur_rel:
                part_of.append((target_id, curr_id))
            if "has_part" in cur_rel:
                has_part.append((curr_id, target_id))
            if "regulates" in cur_rel:
                regulates.append((curr_id, target_id))
        
        #remove redundant edges
        is_a = list(set(is_a))
        part_of = list(set(part_of))
        has_part = list(set(has_part))
        regulates = list(set(regulates))
        
        return (is_a, part_of, has_part, regulates)
