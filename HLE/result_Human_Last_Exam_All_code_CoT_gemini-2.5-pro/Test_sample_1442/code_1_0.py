import networkx as nx
from itertools import combinations

def count_k_matchings(graph, k):
    """
    Counts the number of k-matchings in a given graph.
    A k-matching is a set of k edges with no shared vertices.
    """
    count = 0
    # Iterate through all combinations of k edges from the graph
    for edges in combinations(graph.edges(), k):
        # A set to store all vertices covered by the chosen edges
        vertices = set()
        is_matching = True
        for edge in edges:
            # If a vertex is already in the set, the edges are not disjoint
            if edge[0] in vertices or edge[1] in vertices:
                is_matching = False
                break
            vertices.add(edge[0])
            vertices.add(edge[1])
        
        if is_matching:
            count += 1
    return count

# --- Graph 1: Two disjoint K_3,3 graphs ---
# K_3,3 is a 3-regular bipartite graph on 6 vertices.
# Two disjoint copies give a 3-regular bipartite graph on 12 vertices.
k33_1 = nx.complete_bipartite_graph(3, 3)
k33_2 = nx.complete_bipartite_graph(3, 3)
G1 = nx.disjoint_union(k33_1, k33_2)
G1.name = "Two disjoint K_3,3 graphs"

# --- Graph 2: The hexagonal prism graph ---
# This is another 3-regular bipartite graph on 12 vertices.
G2 = nx.hexagonal_prism_graph()
G2.name = "Hexagonal prism graph"

# The number of edges to choose from (must be the same for both graphs)
# For a d-regular graph on n vertices, |E| = n*d/2
n = 12
d = 3
num_edges = (n * d) // 2

# We are looking for 3-matchings
k = 3

# Calculate the number of 3-matchings for each graph
num_matchings_g1 = count_k_matchings(G1, k)
num_matchings_g2 = count_k_matchings(G2, k)

# Print the results
print(f"Let G1 and G2 be two bipartite, {d}-regular graphs on {n} vertices.")
print(f"G1 = {G1.name}")
print(f"G2 = {G2.name}")
print(f"Both graphs have {num_edges} edges.")
print("-" * 30)
print("The number of 3-matchings is calculated for each graph.")
print(f"Number of 3-matchings in G1: {num_matchings_g1}")
print(f"Number of 3-matchings in G2: {num_matchings_g2}")
print("-" * 30)
if num_matchings_g1 != num_matchings_g2:
    print("The number of 3-matchings is NOT necessarily the same.")
else:
    print("The number of 3-matchings is the same for this example.")
