# -*- coding: utf-8 -*-

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv
import networkx as nx
import matplotlib.pyplot as plt


# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #




# --------------- END OF AUXILIARY FUNCTIONS ------------------ #

def degree_distribution(G: nx.DiGraph, degree_type: str, \
                        node_type: str, bestN: int) -> tuple:
    '''Compute the in-, out- or general degree distribution.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - param: degree_type : str
        Type of degree to compute ('in', 'out', 'general')
    - param: node_type : str
        Type of nodes on which we report degrees ('TF', 'TG', 'all')

    - return: list
        list with degree frequency distribution
    - return: dict
        dictionary with node names as keys and their degrees as values
        for the top N nodes in the degree distribution
    '''

    deg_freq, bestN_degrees = [], {}

    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #

    if G.number_of_nodes()==0:
        return (deg_freq,bestN_degrees)
    nodes_seleccionats=[]
    for node in G.nodes():
        is_selected = False
        
        if node_type == 'all':
            is_selected = True
        elif node_type == 'TF':
            if G.out_degree(node) > 0:
                is_selected = True
        elif node_type == 'TG':
            if G.in_degree(node) > 0:
                is_selected = True
        
        if is_selected:
            nodes_seleccionats.append(node)
    dic_nom = {}
    llista_graus=[]
    for node in nodes_seleccionats:
        if degree_type == 'in':
            dic_nom[node] =G.in_degree(node)
            llista_graus.append((G.in_degree(node)))
            
        elif degree_type == 'out':
            llista_graus.append((G.out_degree(node)))
            dic_nom[node] =G.out_degree(node)
        elif degree_type == 'general':
            llista_graus.append((G.degree(node)))
            dic_nom[node] =G.degree(node)
        
    max_degree = max(llista_graus)
    
    deg_freq = [0] * (max_degree + 1)
    
    for d in llista_graus:
        deg_freq[d] += 1

    sorted_nodes = sorted(dic_nom.items(), key=lambda x: x[1], reverse=True) 
    top_nodes_list = sorted_nodes[:bestN]
    bestN_degrees = dict(top_nodes_list)
    # ----------------- END OF FUNCTION --------------------- #

    return deg_freq, bestN_degrees

def largest_CC_graph(G: nx.Graph) -> nx.Graph:
    '''Generate the largest connected component graph.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - return: Networkx graph
        The graph corresponding to the largest connected component in G
    '''

    H = None
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    if G.is_directed():
        component=max(list(nx.weakly_connected_components(G)))       
    else:
        component=max(list(nx.connected_components(G)))

    if G.number_of_nodes==0:
        return nx.Graph()
    
    H=G.subgraph(component)

    # ----------------- END OF FUNCTION --------------------- #

    return H

if __name__ == "__main__":
    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #
    import time
    start_time = time.time()

    G1=nx.read_graphml('Ecoli_TRN.graphml')
    G2=nx.read_graphml('Ecoli_operon_TRN.graphml')
    #PREGUNTA A
    ordenG1 = G1.order()
    tamaño_G1 = G1.size()
    dis_in_G1 = degree_distribution(G1,'in','all',10)
    dis_out_G1 = degree_distribution(G1,'out','all',10)
    dis_gen_G1 = degree_distribution(G1,'general','all',10)
    component1 = largest_CC_graph(G1)

    

    ordenG2 = G2.order()
    tamaño_G2 = G2.size()
    dis_in_G2 = degree_distribution(G2,'in','all',10)
    dis_out_G2 = degree_distribution(G2,'out','all',10)
    dis_gen_G2 = degree_distribution(G2,'general','all',10)
    component2 = largest_CC_graph(G2)
    print(f" Grafo 1:\n Ordre: {ordenG1}\n Tamaño: {tamaño_G1}\n Distribuició indegree: {dis_in_G1}\n Distribuició outdegree: {dis_out_G1}\n Distribuicio general: {dis_gen_G1}\n Component més llarga: {component1} ")
    print(f" Grafo 1:\n Ordre: {ordenG2}\n Tamaño: {tamaño_G2}\n Distribuició indegree: {dis_in_G2}\n Distribuició outdegree: {dis_out_G2}\n Distribuicio general: {dis_gen_G2}\n Component més llarga: {component2} ")  

    #PREGUNTA B
    dif = []
    for a,b in zip(dis_in_G2[0],dis_in_G1[0]):
        numero = a-b
        dif.append(numero)
    plt.hist(dif, bins=range(max(dif)+2), edgecolor="black")
    plt.xlabel("Grau")
    plt.ylabel("Nombre de nodes")
    plt.title("Histograma de distribució de graus")
    plt.show()
    


    #PREGUNTA C
    # TOP 10 NODES
    indeg_dist_G1, top_in_G1 = degree_distribution(G1, 'in', 'all', 10)
    outdeg_dist_G1, top_out_G1 = degree_distribution(G1, 'out', 'all', 10)

    indeg_dist_G2, top_in_G2 = degree_distribution(G2, 'in', 'all', 10)
    outdeg_dist_G2, top_out_G2 = degree_distribution(G2, 'out', 'all', 10)

    print("\n--- Red SIN operones ---")
    print("Top 10 nodos con mayor IN-degree:")
    for node, val in top_in_G1.items():
        print(f"  {node}  (in-degree = {val}, tipo = {G1.nodes[node].get('ntype')}, gene = {G1.nodes[node].get('name')})")

    print("\nTop 10 nodos con mayor OUT-degree:")
    for node, val in top_out_G1.items():
        print(f"  {node}  (out-degree = {val}, tipo = {G1.nodes[node].get('ntype')}, gene = {G1.nodes[node].get('name')})")

    print("\n--- Red CON operones ---")
    print("Top 10 nodos con mayor IN-degree:")
    for node, val in top_in_G2.items():
        print(f"  {node}  (in-degree = {val}, tipo = {G2.nodes[node].get('ntype')}, gene = {G2.nodes[node].get('name')})")

    print("\nTop 10 nodos con mayor OUT-degree:")
    for node, val in top_out_G2.items():




     
             print(f"  {node}  (out-degree = {val}, tipo = {G2.nodes[node].get('ntype')}, gene = {G2.nodes[node].get('name')})")

    