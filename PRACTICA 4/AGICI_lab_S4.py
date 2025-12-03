# -*- coding: utf-8 -*-
import random
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations


# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #



# --------------- END OF AUXILIARY FUNCTIONS ------------------ #

def largest_CC_graph(G: nx.Graph) -> nx.Graph:
    '''Generate the largest connected component graph.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - return: Networkx graph
        The graph corresponding to the largest connected component in G
    '''

    
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    if G.is_directed():
        component=max(list(nx.weakly_connected_components(G)))        
    else:
        component=max(list(nx.connected_components(G)))

    if G.number_of_nodes==0:
        return nx.Graph(G)
    
    subgraf=G.subgraph(component)
    return subgraf




    # ----------------- END OF FUNCTION --------------------- #

    


def average_distance(G: nx.Graph, iterations: int) -> float:
    '''Estimate the average distance in the input graph.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - param: iterations : int
        Number of iterations to perform
    - return: float
        The estimated average distance in G
    '''

    avg_distance = None
    i=0
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    if G.is_directed():
        G=G.to_undirected()
    
    while i < iterations:
        nodo1=random.choice(list(G.nodes()))
        nodo2=random.choice(list(G.nodes()))
        try:
            avg_distance = nx.shortest_path_length(G, nodo1, nodo2)
        except nx.NetworkXNoPath:
            continue
        
        i += 1

    return avg_distance/iterations





    # ----------------- END OF FUNCTION --------------------- #



def deletion_impact(G: nx.Graph, node_list: list, \
                    grouping_size: int, iterations: int) -> dict:
    '''Assess the impact of node deletions on the graph average distance.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - param: node_list : list
        List of nodes to delete from the network
    - param: grouping_size : list
        The size of the groupings among nodes in the list
    - param: iterations : int
        Number of iterations to perform for average distance
    - return: dict
        Dictionary with grouping node names tuples as keys and differential
        average distance as values.
    '''
    
    del_impact = {}
# ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #

    base_distance = average_distance(G,iterations)

    grupos = combinations(node_list, grouping_size)

    for grupo in grupos:
        G_copy = G.copy()
        G_copy.remove_nodes_from(grupo)

        new_distance = average_distance(G_copy, iterations)

        del_impact[tuple(grupo)] = base_distance - new_distance

    return del_impact




# ----------------- END OF FUNCTION --------------------- #

if __name__ == "__main__":

    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #
    import time
    start_time = time.time()

    G1=nx.read_graphml('Ecoli_TRN.graphml')
    G2=nx.read_graphml('Ecoli_operon_TRN.graphml')
    # APARTAT A
    component=largest_CC_graph(G1)
    avg_dist1=average_distance(G1,10000)
    avg_dist2=nx.average_shortest_path_length(G1)
    print(f"--- DISTANCIA MEDIA 1 (FUNCIÓN)---\n {avg_dist1}")
    print(f"--- DISTANCIA MEDIA 2 (NETWORKX)---\n {avg_dist2}")
    
    
    # APARTAT B

    TF_nodes = [n for n in component.nodes() if component.nodes[n].get('ntype')=='TF']

    impact_individual = deletion_impact(component,TF_nodes, 1, 10000)

    top10 = sorted(impact_individual.items(), key=lambda x: x[1], reverse=True)[:10]
    print(f'\n Top 10 més TF mes impactants')

    for nodes, diff in top10:
        node = nodes[0]
        print(f"{node} - gene: {component.nodes[node].get('name')}, distancia media: {diff:.2f}")
    
    
    # APARTAT C
    top30 = [ tf_tuple[0] for tf_tuple, val in 
          sorted(impact_individual.items(), key=lambda x: x[1], reverse=True)[:30] ]

    impact_pairs = deletion_impact(component, top30, grouping_size=2, iterations=3000)

    top5_pairs = sorted(impact_pairs.items(), key=lambda x: x[1], reverse=True)[:5]

    print("\nTop 5 parelles de TF més impact:")
    for pair, impact_val in top5_pairs:
        print(f"  {pair} --> impact = {impact_val:.4f}")

    print("--- %s seconds ---" % (time.time() - start_time))
    # ------------------- END OF MAIN ------------------------ #




#fer algo de centralitat amb betwenees