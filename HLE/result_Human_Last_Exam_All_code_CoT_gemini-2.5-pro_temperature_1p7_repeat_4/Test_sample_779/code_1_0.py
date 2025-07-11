import networkx as nx

def is_partition_connected(G, partition):
    """Checks if a partition is a 'connected partition' of graph G."""
    for block in partition:
        # A block of size 1 is trivially connected.
        if len(block) > 1:
            # For blocks with more than one vertex, check if the induced subgraph is connected.
            subgraph = G.subgraph(block)
            if not nx.is_connected(subgraph):
                return False
    return True

def is_refinement(p1, p2):
    """Checks if partition p1 is a refinement of partition p2."""
    p1_sets = [set(b) for b in p1]
    p2_sets = [set(b) for b in p2]
    for b1 in p1_sets:
        # Each block in p1 must be a subset of some block in p2
        if not any(b1.issubset(b2) for b2 in p2_sets):
            return False
    return True

# --- Main analysis ---

# 1. Define a simple graph G: a path on 4 vertices
n = 4
G = nx.path_graph(n)
print(f"Graph G has vertices {G.nodes()} and edges {G.edges()}.")

# 2. Define two partitions of the vertices {0, 1, 2, 3}
# Note: networkx uses 0-based indexing by default for path_graph(n)
sigma1 = [{0, 1}, {2}, {3}]
sigma2 = [{0}, {1}, {2, 3}]
print(f"\nConsider partition σ1 = {sigma1}")
print(f"Consider partition σ2 = {sigma2}")

# 3. Verify both partitions are in P(G, n)
is_s1_in_P = is_partition_connected(G, sigma1)
is_s2_in_P = is_partition_connected(G, sigma2)

print(f"\nIs σ1 a connected partition of G? {is_s1_in_P}")
print(f"Is σ2 a connected partition of G? {is_s2_in_P}")

if not (is_s1_in_P and is_s2_in_P):
    print("The chosen partitions are not valid for this example.")
else:
    # 4. Check for comparability to test if P is a total order
    s1_le_s2 = is_refinement(sigma1, sigma2)
    s2_le_s1 = is_refinement(sigma2, sigma1)

    print(f"\nIs σ1 a refinement of σ2 (σ1 <= σ2)? {s1_le_s2}")
    print(f"Is σ2 a refinement of σ1 (σ2 <= σ1)? {s2_le_s1}")

    if not s1_le_s2 and not s2_le_s1:
        print("\nConclusion: σ1 and σ2 are not comparable.")
        print("Therefore, the poset P(G, n) is not a total order.")
        print("This rules out option A.")
    else:
        print("\nThis pair is comparable, try a different pair or graph for a counterexample.")
