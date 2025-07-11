import networkx as nx
import numpy as np

def calculate_query_answers(graph):
    """
    Calculates the number of answers for the query phi(x,y) = exists z. E(x,z) and E(z,y).
    This is equivalent to the number of non-zero elements in the squared adjacency matrix.
    """
    # Get the adjacency matrix of the graph
    # Use nodelist to ensure consistent ordering of vertices
    nodelist = sorted(graph.nodes())
    adj_matrix = nx.to_numpy_array(graph, nodelist=nodelist)
    
    # Square the adjacency matrix to find paths of length 2
    adj_matrix_sq = np.linalg.matrix_power(adj_matrix, 2)
    
    # The number of answers is the number of pairs (u,v) with a path of length 2.
    # This is the number of non-zero entries in the squared adjacency matrix.
    num_answers = np.count_nonzero(adj_matrix_sq)
    
    return adj_matrix_sq, num_answers

# --- Graph 1: Path graph on 6 vertices (P_6) ---
G1 = nx.path_graph(6)
G1_adj_sq, G1_answers = calculate_query_answers(G1)

print("Graph G1: Path graph on 6 vertices (P_6)")
print("Degree sequence:", [d for n, d in G1.degree()])
print("Squared Adjacency Matrix (A_G1^2):")
print(G1_adj_sq.astype(int))
print(f"Number of answers for phi in G1: {G1_answers}\n")


# --- Graph 2: Disjoint union of a 4-cycle and an edge (C_4 U K_2) ---
G2 = nx.disjoint_union(nx.cycle_graph(4), nx.path_graph(2))
# Relabel nodes for clarity if desired, but not necessary for the calculation.
# For example: mapping = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f'}
# G2 = nx.relabel_nodes(G2, mapping)
G2_adj_sq, G2_answers = calculate_query_answers(G2)

print("Graph G2: Disjoint union of C_4 and K_2")
print("Degree sequence:", sorted([d for n, d in G2.degree()]))
print("Squared Adjacency Matrix (A_G2^2):")
print(G2_adj_sq.astype(int))
print(f"Number of answers for phi in G2: {G2_answers}\n")

# --- Conclusion ---
if G1_answers != G2_answers:
    print("The number of answers is different for G1 and G2.")
    print(f"Therefore, it is possible that G1 and G2 have different numbers of answers for the query phi.")
else:
    print("The number of answers is the same for G1 and G2 for this query.")
