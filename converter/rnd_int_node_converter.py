"""
To run this code:
python3 converter_updated.py -i "INPUT FILE PATH"

Il criterio sceglie un nodo interno "casualmente", cioè l'ultimo nodo interno che incontra è quello che viene tenuto
(in realtà è casuale se gli archi degli alberi sono inseriti a caso e senza seguire uno schema logico preciso)
TODO
- in questo converter come criterio di scelta di nodo significativo scelgo i nodi interni,
quindi nodi radice e nodi padre sono meno significativi

"""

import numpy as np
import pandas as pd
from itertools import combinations
from collections import defaultdict
import argparse

# Create parser
parser = argparse.ArgumentParser()

# adds arguments
parser.add_argument("-i", help="Input file path")

args = parser.parse_args()
input_file = args.i
output_file = (args.i).replace(".txt", "")  + "-converted2.txt"

df = pd.read_csv('data/table_gene_path.csv', delimiter=';', header=None, encoding='utf-8')

data_dict = {}

# Iterazione sul dataframe per elaborare ogni riga
for index, row in df.iterrows():
    gene = row[0]
    pathways = [pathway.strip() for pathway in row[1].split(';')]
    
    # Salva i dati nel dizionario
    data_dict[gene] = pathways


sorted_keys = sorted(data_dict.keys(), key=len, reverse=True)

# custom names to pathways 
df_custom = pd.read_csv('data/pathways.csv', delimiter=';', header=None, encoding='utf-8')
replacement_dict = dict(zip(df_custom[0], df_custom[1]))


for key, values in data_dict.items():
    # Sostituisci i valori nella lista in base al dizionario di sostituzione
    data_dict[key] = [replacement_dict.get(value, value) for value in values]


#print(data_dict)



# types of edges
dir_edge, diff_edge, un_edge = "->-", r"-/-", "-?-" 
# dict to store score of each node
score_dict = {}

output_trees = ""



# update score (amp handled)
"""  
- root has descentant set to 1
- leafs have ancestor set to 1
- internal nodes have both ancestor and desccentant set to one
TODO capire cosa succede se ho due branch e ci sono due nodi interni con lo stesso pathway
"""
def update_score(node_dict, edge):
    # if dir_edge update the scores
    if dir_edge in edge:
        nodes = edge.split(dir_edge)
        first_node, second_node = nodes[0], nodes[1]

        # first_node of the edge      
        if first_node not in node_dict:
            if 'amp' in first_node:
                node_dict[first_node.replace("amp", "")] = {'pathways': data_dict[first_node.replace("amp", "")], 'ancestor': 0, 'descentant': 1, 'amp': 1}
            else: 
                node_dict[first_node.replace("amp", "")] = {'pathways': data_dict[first_node.replace("amp", "")], 'ancestor': 0, 'descentant': 1, 'amp': 0}
        else:
           node_dict[first_node]['descentant'] = 1
        
        # second_node of the edge
        if second_node not in node_dict:
            if 'amp' in second_node:
                node_dict[second_node.replace("amp", "")] = {'pathways': data_dict[second_node.replace("amp", "")], 'ancestor': 1, 'descentant': 0, 'amp': 1}
            else:
                node_dict[second_node.replace("amp", "")] = {'pathways': data_dict[second_node.replace("amp", "")], 'ancestor': 1, 'descentant': 0, 'amp': 0}
        else:
            node_dict[second_node]['ancestor'] = 1
    
    # if diff_edge just add to the list if not already in it 
    elif diff_edge in edge:
        nodes = edge.split(diff_edge)
        first_node, second_node = nodes[0], nodes[1]

        if first_node not in node_dict:
            if 'amp' in first_node:
                node_dict[first_node.replace("amp", "")] = {'pathways': data_dict[first_node.replace("amp", "")], 'ancestor': 0, 'descentant': 0, 'amp': 1}
            else: 
                node_dict[first_node.replace("amp", "")] = {'pathways': data_dict[first_node.replace("amp", "")], 'ancestor': 0, 'descentant': 0, 'amp': 0}
        if second_node not in node_dict:
            if 'amp' in second_node:
                node_dict[second_node.replace("amp", "")] = {'pathways': data_dict[second_node.replace("amp", "")], 'ancestor': 0, 'descentant': 0, 'amp': 1}
            else:
                node_dict[second_node.replace("amp", "")] = {'pathways': data_dict[second_node.replace("amp", "")], 'ancestor': 0, 'descentant': 0, 'amp': 0}
        
    elif un_edge:
        nodes = edge.split(un_edge)
        first_node, second_node = nodes[0], nodes[1]

        if first_node not in node_dict:
            if 'amp' in first_node:
                node_dict[first_node.replace("amp", "")] = {'pathways': data_dict[first_node.replace("amp", "")], 'ancestor': 0, 'descentant': 0, 'amp': 1}
            else: 
                node_dict[first_node.replace("amp", "")] = {'pathways': data_dict[first_node.replace("amp", "")], 'ancestor': 0, 'descentant':0, 'amp': 0}
        if second_node not in node_dict:
            if 'amp' in second_node:
                node_dict[second_node.replace("amp", "")] = {'pathways': data_dict[second_node.replace("amp", "")], 'ancestor': 0, 'descentant': 0, 'amp': 1}
            else:
                node_dict[second_node.replace("amp", "")] = {'pathways': data_dict[second_node.replace("amp", "")], 'ancestor': 0, 'descentant': 0, 'amp': 0}
        

def create_edge(converter, first_node, second_node, edge, temp):
    first_node = first_node.replace(first_node, '; '.join(converter[first_node]))
    second_node = second_node.replace(second_node, '; '.join(converter[second_node]))
    # if there are more pathways in one node
    first_node = first_node.split("; ")
    second_node = second_node.split("; ")

    # metto -?- per tutti i pathway nello stesso nodo
    combinazioni = list(combinations(first_node, 2))
    
    for combo in combinazioni:
        if f"{combo[0]}{un_edge}{combo[1]}" not in temp:
            temp += f"{combo[0]}{un_edge}{combo[1]} "

    combinazioni = list(combinations(second_node, 2))
    for combo in combinazioni:
        if f"{combo[0]}{un_edge}{combo[1]}" not in temp:
            temp += f"{combo[0]}{un_edge}{combo[1]} "
    
    # make combinations only if both nodes have values
    if first_node != [''] and second_node != ['']:
        # metto ->- per i pathway contigui
        X, Y = np.meshgrid(first_node, second_node)
        matrix = np.stack((X,Y), axis=-1)
        # Iterate through each element in the matrix
        for row in matrix:
            for pair in row:
                # Print the pair in the desired format
                if f"{pair[0]}{edge}{pair[1]}" not in temp:
                    temp += f"{pair[0]}{edge}{pair[1]} " 

    return temp

def pathway_importance_criterio(node_dict):
    pathway_info = defaultdict(lambda: {'importance': 0,'gene': ''})
    for gene, info in node_dict.items():
        for pathway in info['pathways']:
            if info['ancestor'] == 1 and info['descentant'] == 1:
                pathway_info[pathway]['importance'] = 1
                pathway_info[pathway]['gene'] = gene
            elif (info['ancestor'] == 0 or info['descentant'] == 0) and pathway not in pathway_info:
                pathway_info[pathway]['importance'] = 0
                pathway_info[pathway]['gene'] = gene
    return pathway_info

 
# three phases of the algorithm
# might 2 + 3 can be a single phase
with open(input_file, 'r') as file:
    for tree in file:

        # update score phase
        node_dict = {}
        tree = tree.strip().split()
        for edge in tree:  
            update_score(node_dict, edge)

        print(node_dict)        
        print("\n")

        # pathways importance criterio 
        # criterio if ancestor=1 and descentant=1 is more imporant than {a=1, d=0} or {a=0, d=1}
        
        pathway_info = pathway_importance_criterio(node_dict)    

        # filtered dict with importance criterio
        #print(node_dict)
        node_dict_filtered = {}
        
        # handled gene with 'amp'
        for gene, info in node_dict.items():
            updated_pathways = []
            if info['amp'] == 1:
                for pathway in info['pathways']:
                    if gene == pathway_info[pathway]['gene']:
                        updated_pathways.append(pathway+'amp')
                node_dict_filtered[gene +'amp'] = updated_pathways
            else: 
                for pathway in info['pathways']:
                    if gene == pathway_info[pathway]['gene']:
                        updated_pathways.append(pathway)
                node_dict_filtered[gene] = updated_pathways
        

        print(node_dict_filtered)
    # substituting phase
    # questa parte la ho già fatta nel codice vecchio ed è l'unica che funziona basta metterla apposto in modo che funzioni in generale
        temp = ""

        for edge in tree:
            # if ->- is the edge 
            if dir_edge in edge:
                nodes = edge.split(dir_edge)
                temp = create_edge(node_dict_filtered, nodes[0], nodes[1], dir_edge, temp) 
            #if -?- is the edge 
            elif un_edge in edge:
                nodes = edge.split(un_edge)
                temp = create_edge(node_dict_filtered, nodes[0], nodes[1], un_edge, temp) 
            elif diff_edge in edge:
                nodes = edge.split(diff_edge)
                temp = create_edge(node_dict_filtered, nodes[0], nodes[1], diff_edge, temp) 
        
        # if tree doesn't exists output has already \n in it
        if temp != "":
            temp += "\n" # end of the current tree
        output_trees += temp # add tree to output file


with open(output_file, "w") as file:
    file.write(output_trees)