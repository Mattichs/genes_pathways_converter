from matplotlib.ticker import MaxNLocator
import numpy as np
import argparse
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with graphs")
args = parser.parse_args()

reading_ = True
seps = ["->-" , "-/-" , "-?-"]
graphs_list = list()

def longest_path_in(G):
    # Ottieni l'ordinamento topologico del grafo
    topological_order = list(nx.topological_sort(G))
    
    # Inizializza le distanze con -inf e il nodo di partenza con 0
    dist = {node: float('-inf') for node in G.nodes()}
    dist[topological_order[0]] = 0
    
    # Calcola le distanze nel'ordinamento topologico
    for node in topological_order:
        for succ in G.successors(node):
            if dist[succ] < dist[node] + 1:  # Assumiamo pesi unitari
                dist[succ] = dist[node] + 1
    
    # Trova il nodo con la distanza massima
    max_dist = max(dist.values())
    
    return max_dist

def print_graph(tree):
    nx.draw(tree,with_labels=True)
    plt.draw()
    plt.show()

def load_tree(line):
    edges_ = line.split(" ")
    G = nx.DiGraph()
    G_tree = nx.DiGraph()
    # first compute the nodes of the tree
    for edge in edges_:
        if "-?-" in edge:
            # these two alterations go in the same node
            nodes = edge.split("-?-")
            G.add_edge(nodes[0] , nodes[1])
    
    # each connected component of G is a node of the tree
    components = nx.weakly_connected_components(G)
    for comp in components:
        list_nodes_comp = [i for i in comp]
        list_nodes_comp.sort()
        comp_str = list_nodes_comp[0]
        for i in range(1,len(list_nodes_comp)):
            comp_str = comp_str+ "-"+ list_nodes_comp[i]
        G_tree.add_node(comp_str)
    for edge in edges_:
        if "-/-" in edge:
            nodes = edge.split("-/-")
            node_1 = nodes[0]
            node_2 = nodes[1]
            for node_ in G_tree.nodes:
                if nodes[0] in node_:
                    node_1 = node_
                if nodes[1] in node_:
                    node_2 = node_
            if node_1 not in G_tree.nodes:
                G_tree.add_node(node_1)
            if node_2 not in G_tree.nodes:
                G_tree.add_node(node_2)

        if "->-" in edge:
            nodes = edge.split("->-")
            node_1 = nodes[0]
            node_2 = nodes[1]
            for node_ in G_tree.nodes:
                if nodes[0] in node_:
                    node_1 = node_
                if nodes[1] in node_:
                    node_2 = node_
            G_tree.add_edge(node_1 , node_2)

    G_tree.add_node("g")
    for node in G_tree.nodes:
        if node != "g":
            G_tree.add_edge("g" , node)

    return G_tree

def compute_avg_nodes(fin_graphs): 
    debug_loading_graph = 0
    num_nodes = dict()
    i = 1
    for line in fin_graphs:
        line = line.replace("\n","")
        graph_ = load_tree(line) # for every line in the file load tree
        if debug_loading_graph == 1:
            print("loaded shown graph")
            print("line",line)
            print_graph(graph_)
        graphs_list.append(graph_)

        num_nodes_tree = len(list(graph_.nodes))

        if num_nodes_tree in num_nodes:
            num_nodes[num_nodes_tree] = num_nodes[num_nodes_tree] + 1
        else:
            num_nodes[num_nodes_tree] = 1

    print("loaded",len(graphs_list),"graphs")

    # average number of nodes 
    avg_numnodes = 0.
    # average depth of the graphs
    avg_depth = 0.
    # total number of the graphs
    tot_graphs = sum(list(num_nodes.values()))
    
    for numnode_ in num_nodes:
        avg_numnodes += (numnode_)*num_nodes[numnode_]/tot_graphs
    print("average num_nodes: ", avg_numnodes)

    return num_nodes, avg_numnodes

def compute_avg_depth(fin_graphs): 
    debug_loading_graph = 0
    num_depth = dict()
    for line in fin_graphs:
        line = line.replace("\n","")
        graph_ = load_tree(line) # for every line in the file load tree
        if debug_loading_graph == 1:
            print("loaded shown graph")
            print("line",line)
            print_graph(graph_)
        graphs_list.append(graph_)

        # find depth of the graph

        depth = longest_path_in(graph_)
        if depth in num_depth: 
            num_depth[depth] = num_depth[depth] + 1
        else:
            num_depth[depth] = 1

    print("loaded",len(graphs_list),"graphs")

    
    # average depth of the graphs
    avg_depth = 0.
    # total number of the graphs
    tot_graphs = sum(list(num_nodes.values()))

    for depth_ in num_depth:
        avg_depth += (depth_)*num_depth[depth_]/tot_graphs
    print("average num_depth: ", avg_depth)

    return num_depth, avg_depth

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for ax in axs:
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

fin_graphs_1 = open(args.i,"r")

num_nodes, avg_nodes = compute_avg_nodes(fin_graphs_1)
# Primo grafico
width = 0.5
axs[0].bar(np.array(list(num_nodes.keys())), num_nodes.values(), width)#, color='blue')
axs[0].set_title(f"Genes\nDistribution of number of nodes (avg={str(avg_nodes)[0:4]})")
axs[0].set_xlabel("Number of nodes")
axs[0].set_ylabel("Number of trees")
axs[0].set_ylim(0, 70) 


fin_graphs_2 = open(f"{(args.i).replace(".txt", "")}-converted1.txt", "r")

num_nodes, avg_nodes = compute_avg_nodes(fin_graphs_2)
# Secondo grafico
width = 0.5
axs[1].bar(np.array(list(num_nodes.keys())), num_nodes.values(), width)#, color='blue')
axs[1].set_title(f"Pathways 1\nDistribution of number of nodes (avg={str(avg_nodes)[0:4]})")
axs[1].set_xlabel("Number of nodes")
axs[1].set_ylabel("Number of trees")
axs[1].set_ylim(0, 70) 



fin_graphs_3 = open(f"{(args.i).replace(".txt", "")}-converted2.txt", "r")
num_nodes, avg_nodes = compute_avg_nodes(fin_graphs_3)
# Terzo grafico
width = 0.5
axs[2].bar(np.array(list(num_nodes.keys())), num_nodes.values(), width)#, color='blue')
axs[2].set_title(f"Pathways 2\nDistribution of number of nodes (avg={str(avg_nodes)[0:4]})")
axs[2].set_xlabel("Number of nodes")
axs[2].set_ylabel("Number of trees")
axs[2].set_ylim(0, 70) 


#plt.show()
plt.savefig('analisi_alberi_nodes.png', bbox_inches='tight')


#fig.clear()
# ANALISI PROFONDITA'

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for ax in axs:
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

fin_graphs_1 = open(args.i,"r")
num_depth, avg_depth  = compute_avg_depth(fin_graphs_1)
# Primo grafico
width = 0.5
axs[0].bar(np.array(list(num_depth.keys())), num_depth.values(), width)#, color='blue')
axs[0].set_title(f"Genes\nDistribution of Depth (avg={str(avg_depth)[0:4]})")
axs[0].set_xlabel("Depth")
axs[0].set_ylabel("Number of trees")
axs[0].set_ylim(0, 70) 


fin_graphs_2 = open(f"{(args.i).replace(".txt", "")}-converted1.txt", "r")
num_depth, avg_depth = compute_avg_depth(fin_graphs_2)
# Secondo grafico
width = 0.5
axs[1].bar(np.array(list(num_depth.keys())), num_depth.values(), width)#, color='blue')
axs[1].set_title(f"Pathways 1\nDistribution of Depth(avg={str(avg_depth)[0:4]})")
axs[1].set_xlabel("Depth")
axs[1].set_ylabel("Number of trees")
axs[1].set_ylim(0, 70) 


fin_graphs_3 = open(f"{(args.i).replace(".txt", "")}-converted2.txt", "r")
num_depth, avg_depth = compute_avg_depth(fin_graphs_3)
# Terzo grafico
width = 0.5
axs[2].bar(np.array(list(num_depth.keys())), num_depth.values(), width)#, color='blue')
axs[2].set_title(f"Pathways 2\nDistribution of Depth (avg={str(avg_depth)[0:4]})")
axs[2].set_xlabel("Depth")
axs[2].set_ylabel("Number of trees")
axs[2].set_ylim(0, 70) 


#plt.show()
plt.savefig('analisi_alberi_depth.png', bbox_inches='tight')




# ANALISTI TOT. PATTERN E MEDIA SUPPORTI
import pandas as pd

fin_final = open(f"{(args.i).replace(".txt", "")}_final.txt", "r")
df = pd.read_csv(fin_final, sep=';')
print(df["traj_supp"])
tot_trees = df.shape[0] 
traj = df["traj_supp"].sum()
data = {
    "Gene":  [tot_trees, traj/tot_trees]
}
indici = ['Nr. Pattern', 'Media Traiettorie']
table = pd.DataFrame(data, index=indici)

fin_final = open(f"{(args.i).replace(".txt", "")}-converted1_final.txt", "r")
df = pd.read_csv(fin_final, sep=';')
tot_trees = df.shape[0] 
traj = df["traj_supp"].sum()
table["Pathways 1"] = [tot_trees, traj/tot_trees]


fin_final = open(f"{(args.i).replace(".txt", "")}-converted2_final.txt", "r")
df = pd.read_csv(fin_final, sep=';')
tot_trees = df.shape[0] 
traj = df["traj_supp"].sum()
table["Pathways 2"] = [tot_trees, traj/tot_trees]
# Creazione del grafico della tabella con matplotlib
plt.figure(figsize=(8, 4))  # Imposta la dimensione della figura
plt.axis('off')  # Nasconde gli assi per ottenere solo la tabella

# Disegna la tabella con pandas
out_table = plt.table(cellText=table.values,
                  colLabels=table.columns,
                  rowLabels=table.index,  # Aggiunge le etichette di riga
                  loc='center')

# Imposta il layout della tabella
out_table.auto_set_font_size(False)
out_table.set_fontsize(12)
out_table.scale(1.5, 1.5)
plt.savefig('tabella.png', bbox_inches='tight')
plt.show()