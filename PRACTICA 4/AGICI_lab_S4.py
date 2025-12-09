# -*- coding: utf-8 -*-
import random
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations


# ------- AUXILIARY FUNCTIONS ------- #



# --------------- END OF AUXILIARY FUNCTIONS ------------------ #


def largest_CC_graph(G: nx.Graph) -> nx.Graph:
    if G.is_directed():
        component = max(list(nx.weakly_connected_components(G)))
    else:
        component = max(list(nx.connected_components(G)))

    if G.number_of_nodes() == 0:
        return nx.Graph(G)

    return G.subgraph(component)


def average_distance(G: nx.Graph, iterations: int) -> float:
    if G.is_directed():
        G = G.to_undirected()

    nodes = list(G.nodes())
    total = 0
    count = 0

    while count < iterations:
        a = random.choice(nodes)
        b = random.choice(nodes)
        if a == b:
            continue
        try:
            d = nx.shortest_path_length(G, a, b)
        except nx.NetworkXNoPath:
            continue

        total += d
        count += 1

    return total / count


def deletion_impact(G: nx.Graph, node_list: list,
                    grouping_size: int, iterations: int) -> dict:

    del_impact = {}
    base = average_distance(G, iterations)

    for group in combinations(node_list, grouping_size):
        Gc = G.copy()
        Gc.remove_nodes_from(group)
        new = average_distance(Gc, iterations)
        del_impact[tuple(group)] = base - new

    return del_impact


# ===================== EXERCICI 3A ===================== #
def impact_key(x):
    return x[1]


if __name__ == "__main__":

    G = nx.read_graphml("Ecoli_TRN.graphml")

    CC = largest_CC_graph(G).to_undirected()

    print("Nodes CC:", CC.number_of_nodes())
    print("Arestes CC:", CC.number_of_edges())

    exact = nx.average_shortest_path_length(CC)
    print("Distància exacta:", exact)

    iters_list = [10, 100, 1000, 5000, 10000]
    for it in iters_list:
        print(it, "iteracions ->", average_distance(CC, it))

    N = 10000
    n = CC.number_of_nodes()
    total_pairs = n * (n - 1) // 2
    print("Parelles exactes:", total_pairs)
    print("Speed-up aprox:", total_pairs / N)


    # ===================== EXERCICI 3B ===================== #

    tf_nodes = [node for node, data in CC.nodes(data=True)
                if data.get("ntype") == "TF"]

    print("TF totals:", len(tf_nodes))

    impact1 = deletion_impact(CC, tf_nodes, 1, N)
    ordered1 = sorted(impact1.items(), key=impact_key, reverse=True)

    print("Top 10 TF eliminats:")
    for group, diff in ordered1[:10]:
        name = CC.nodes[group[0]].get("name", "")
        print(group[0], name, "impacte:", diff)


    # ===================== EXERCICI 3C ===================== #

    top30 = [g[0] for g in ordered1[:30]]
    impact2 = deletion_impact(CC, top30, 2, N)
    ordered2 = sorted(impact2.items(), key=impact_key, reverse=True)

    print("Top 5 parelles:")
    for group, diff in ordered2[:5]:
        a, b = group
        n1 = CC.nodes[a].get("name", "")
        n2 = CC.nodes[b].get("name", "")
        print(a, n1, "-", b, n2, "impacte:", diff)


    # ===================== EXERCICI 3D ===================== #

    best_TF = ordered1[0][0][0]
    best_name = CC.nodes[best_TF].get("name", "")

    print("\n--- Raonament exercici 3d ---")
    print("El TF que triaria per deshabilitar parcialment la cèl·lula és:", best_TF, best_name)
    print("És el que té l'impacte més alt en augmentar la distància mitjana de la xarxa,")
    print("fet que indica que és estructuralment essencial per mantenir la regulació global.")
    print("Un antibiòtic que ataqui aquest TF podria afectar molts gens alhora,")
    print("reduint la capacitat d'adaptació de la cèl·lula i dificultant l'aparició de resistència.")
