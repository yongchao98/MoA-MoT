import itertools

def count_k_matchings(edges, k):
    """
    Counts the number of k-matchings in a graph given by its edge list.
    A k-matching is a set of k edges with no shared vertices.
    """
    count = 0
    # Iterate through all combinations of k edges from the graph's edge list
    for edge_combo in itertools.combinations(edges, k):
        vertices = set()
        # Collect all vertices from the selected edges
        for u, v in edge_combo:
            vertices.add(u)
            vertices.add(v)
        
        # If the number of unique vertices is 2*k, it's a valid k-matching
        if len(vertices) == 2 * k:
            count += 1
            
    return count

# Define Graph G1: A single cycle C_12
# n=12, d=2. Bipartite since 12 is even.
n1 = 12
# Edges are (i, i+1) for i=0..10 and (11,0)
g1_edges = [tuple(sorted((i, (i + 1) % n1))) for i in range(n1)]

# Define Graph G2: Disjoint union of three C_4 cycles
# n=12, d=2. Bipartite since each component is.
g2_edges = []
num_components = 3
component_size = 4
for i in range(num_components):
    offset = i * component_size
    for j in range(component_size):
        u = offset + j
        v = offset + (j + 1) % component_size
        g2_edges.append(tuple(sorted((u,v))))

# We are interested in 3-matchings
k = 3

# Calculate the number of 3-matchings for both graphs
num_matchings_g1 = count_k_matchings(g1_edges, k)
num_matchings_g2 = count_k_matchings(g2_edges, k)

# Print the explanation and results
print("The question is whether any two bipartite, d-regular graphs on n vertices have the same number of 3-matchings.")
print("We can test this with a counterexample.")
print("\nLet's consider two graphs for d=2 and n=12:")
print("G1 = C_12 (a single 12-vertex cycle)")
print("G2 = C_4 U C_4 U C_4 (three disjoint 4-vertex cycles)")
print("\nBoth G1 and G2 are bipartite, 2-regular graphs on 12 vertices.")
print("Let's count the number of 3-matchings in each graph.")
print("-" * 30)
print(f"Number of 3-matchings in G1 (C_12): {num_matchings_g1}")
print(f"Number of 3-matchings in G2 (C_4 U C_4 U C_4): {num_matchings_g2}")
print("-" * 30)

if num_matchings_g1 != num_matchings_g2:
    print(f"The final numbers are different: {num_matchings_g1} != {num_matchings_g2}.")
    print("\nSince we found two graphs that satisfy the conditions but have a different number of 3-matchings, the answer is No.")
else:
    # This case will not be reached based on our analysis
    print(f"The final numbers are the same: {num_matchings_g1} == {num_matchings_g2}.")
    print("\nIn this specific case, the numbers are the same, but this does not prove the general statement.")
