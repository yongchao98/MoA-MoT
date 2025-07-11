import itertools

def find_k_cliques(graph, k):
    """
    Finds all k-cliques in a given graph.
    A k-clique is a set of k vertices where every two vertices are connected.
    
    Args:
        graph (dict): An adjacency list representation of the graph.
        k (int): The size of the clique to find.
        
    Returns:
        list: A list of lists, where each inner list represents a k-clique.
    """
    # This algorithm demonstrates the #k-Clique problem by brute-force enumeration.
    # Its complexity is roughly O(n^k * k^2), which is not FPT.
    nodes = list(graph.keys())
    if k > len(nodes):
        return []
    
    cliques = []
    # Step 1: Generate all combinations of k vertices from the graph
    for vertex_combo in itertools.combinations(nodes, k):
        is_clique = True
        # Step 2: For each combination, check if all pairs of vertices are connected
        for u, v in itertools.combinations(vertex_combo, 2):
            # Adjacency check using the graph representation
            if v not in graph.get(u, []):
                is_clique = False
                break
        
        # Step 3: If all pairs are connected, it's a clique
        if is_clique:
            cliques.append(list(vertex_combo))
            
    return cliques

# --- Example for the problem #Sub_G(H) where H is the class of cliques ---
# Input: A graph G and an integer k (since H is K_k)
# G is a sample graph from a somewhere dense class (the class of all graphs)
G_adj = {
    0: [1, 2, 3, 4],
    1: [0, 2, 3, 4],
    2: [0, 1, 3],
    3: [0, 1, 2],
    4: [0, 1]
}

# Parameter: |H| = k. Let's find all 3-cliques (triangles).
k = 3

# Output: The number of subgraphs of G isomorphic to K_k
found_cliques = find_k_cliques(G_adj, k)
count = len(found_cliques)

print(f"This script demonstrates the problem from option B: counting k-cliques.")
print(f"Input graph G has vertices: {list(G_adj.keys())}")
print(f"Pattern H is a clique of size k = {k}.")
print("-" * 20)
print(f"Found {count} cliques of size {k}:")
for clique in found_cliques:
    print(f"- {clique}")

# The prompt requires printing numbers in a final equation.
# We can represent the total count as a sum.
equation_parts = ["1"] * count
equation_str = " + ".join(equation_parts) if count > 0 else "0"

print("\nThe final count is derived from summing 1 for each clique found:")
# Final equation: 1 + 1 + ... + 1 = total_count
if count > 0:
    print(f"{equation_str} = {count}")
else:
    print(f"Count = 0")
