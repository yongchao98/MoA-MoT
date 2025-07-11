import itertools

def count_subgraphs_isomorphic_to_clique(G, k):
    """
    This function counts the number of k-cliques in a graph G.
    This is the problem #Sub_G(H) where H is a k-clique.
    This brute-force implementation illustrates the combinatorial nature of the problem.
    Its runtime complexity is roughly O(n^k), where n is the number of vertices,
    which is not fixed-parameter tractable.

    Args:
        G (dict): An adjacency list representation of the graph.
                  e.g., {0: [1, 2], 1: [0], 2: [0]}
        k (int): The size of the clique to count.

    Returns:
        int: The number of k-cliques found in G.
    """
    nodes = list(G.keys())
    n = len(nodes)
    if k > n:
        return 0
    
    clique_count = 0
    
    # Iterate over all combinations of k vertices
    for subset_indices in itertools.combinations(range(n), k):
        # Extract the node labels for the current subset
        subset_nodes = [nodes[i] for i in subset_indices]
        is_clique = True
        
        # Check all pairs in the subset for an edge
        for i in range(k):
            for j in range(i + 1, k):
                u, v = subset_nodes[i], subset_nodes[j]
                # If any pair is not connected, it's not a clique
                if u not in G[v]:
                    is_clique = False
                    break
            if not is_clique:
                break
        
        if is_clique:
            clique_count += 1
            
    return clique_count

# --- Example ---
# Statement B asserts that for H = cliques, the problem is #W[1]-complete.
# Let's illustrate the counting with a concrete example.
# Let G be a complete graph K_5 (a 5-clique).
# Let H be a K_3 (a triangle).
# The number of subgraphs of K_5 isomorphic to K_3 is "5 choose 3".
# Equation: C(5, 3) = 10

# Graph G = K_5
K5_graph = {
    0: [1, 2, 3, 4],
    1: [0, 2, 3, 4],
    2: [0, 1, 3, 4],
    3: [0, 1, 2, 4],
    4: [0, 1, 2, 3]
}
k_size = 3

# Calculate the number using the function
number_of_triangles = count_subgraphs_isomorphic_to_clique(K5_graph, k_size)

print("Illustrating the counting problem for statement B:")
print(f"Let host graph G be a K_5 (clique on 5 vertices).")
print(f"Let pattern graph H be a K_3 (clique on 3 vertices).")
print(f"The number of subgraphs of G isomorphic to H is: {number_of_triangles}")
print("\nThis corresponds to the equation for combinations C(n, k) = result, or C(5, 3) = 10.")
print("The numbers in this equation are:")
print(5)
print(3)
print(10)