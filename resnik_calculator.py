import networkx as nx 
import numpy as np
import numpy as np
import pandas as pd
import os
import sys


'''
=== Need ===
go: Gene Ontology table with three columns: GO_term, term_size, genes (in Kramer format)
genelist: list of genes to calculate resnik similarity for
N: total number of genes
'''
class resnik_calculator:
    def __init__(self, go_url, gene_list):
        self.go = pd.read_table(go_url)
        self.go.columns = ['term', 'tsize', 'genes']
        self.N = self.go['tsize'].max()
        self.genelist = gene_list

    def get_resnik(self):
        # Sort GO by term size
        self.go.sort_values('tsize', inplace=True)
        # Calculate resnik score
        resnik = pd.DataFrame(-1, index=self.genelist, columns=self.genelist, dtype=float)
        geneset = set(self.genelist)
        count = 0 # for process tracking

        for idx, row in self.go.iterrows():
            count += 1
            if count % 10 == 0:
                print('... processed {}%'.format(count*100/self.go.shape[0]))
            rsim = -np.log10(row['tsize'] / self.N) # resnik similarity
            genes = set(row['genes'].split(',')[:-1]).intersection(geneset)
            if len(genes) > 1:
                genes = list(genes)
                for i in range(len(genes)-1):
                    ga = genes[i]
                    for j in range(i+1, len(genes)):
                        gb = genes[j]
                        # fill in the resnik similarity only if gene pair is unseen before
                        if resnik.at[ga, gb] == -1:
                            resnik.at[ga, gb] = rsim
                            resnik.at[gb, ga] = rsim

        # Assign genes never appeared in the same term similarity of 0
        resnik.replace(-1, 0, inplace=True)
        # Fill in the diagonal with maximum resnik value
        max_resnik = resnik.max().max()
        np.fill_diagonal(resnik.values, max_resnik)
        # Scale resnik similarity into [0,1]
        print('Max raw resnik is {}'.format(resnik.max().max()))
        print('Min raw resnik is {}'.format(resnik.min().min()))
        resnik_scaled = resnik / max_resnik
        return (resnik, resnik_scaled)
  

rc = resnik_calculator('./human_go_cc.no_hpa.symbol.termStats', ['TP53', 'CCNB1', 'DHX3', 'GSDMD', 'FBXL12']) 
r, r_s = rc.get_resnik()
print(r)
print(r_s)