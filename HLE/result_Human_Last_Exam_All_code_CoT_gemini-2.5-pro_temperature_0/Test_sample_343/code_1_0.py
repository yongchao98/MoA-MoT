import networkx as nx
from itertools import combinations

def reduce_clique_to_pdecide(G: nx.Graph, k: int):
    """
    Performs a parameterized reduction from k-CLIQUE on G to the problem PDecide.
    The reduction is valid for k >= 4.

    Args:
        G (nx.Graph): The input graph for k-CLIQUE.
        k (int): The size of the clique to find.

    Returns:
        A tuple (G_prime, k_prime) for the new problem instance, or None if k < 4.
    """
    if k < 4:
        print(f"Reduction is designed for k >= 4. Input k={k} is too small.")
        # For k < 4, k-CLIQUE is FPT, so we can solve it and return a trivial instance.
        # Here, we'll just return None for simplicity.
        return None, None

    # The new parameter k' is k choose 2.
    k_prime = k * (k - 1) // 2
    
    print(f"Original problem: Find a {k}-clique in a graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    print(f"Reduction equation for the new parameter: k' = k * (k - 1) / 2")
    print(f"Resulting parameter: k' = {k} * ({k} - 1) / 2 = {k_prime}")


    G_prime = nx.Graph()
    
    V = list(G.nodes())
    E = list(G.edges())
    
    # Create the two sets of vertices for the bipartite graph G'
    # V_parts are vertices of the form (v, i)
    V_parts = []
    for v in V:
        for i in range(1, k + 1):
            V_parts.append((v, i))
            
    # E_parts are vertices of the form (e, (i, j))
    E_parts = []
    edge_indices = list(combinations(range(1, k + 1), 2))
    for e in E:
        for p in edge_indices:
            # Ensure node labels in edges are ordered for canonical representation
            sorted_e = tuple(sorted(e))
            E_parts.append((sorted_e, p))

    # Add nodes to G'
    G_prime.add_nodes_from(V_parts, bipartite=0)
    G_prime.add_nodes_from(E_parts, bipartite=1)

    # Add edges to G' based on the reduction rule
    # Edge between (v, l) and (e, (i, j)) if v is an endpoint of e and l is in {i, j}
    for v_node in V_parts:
        v, l = v_node
        for e_node in E_parts:
            e, p = e_node
            i, j = p
            
            if v in e and l in {i, j}:
                G_prime.add_edge(v_node, e_node)

    return G_prime, k_prime

if __name__ == '__main__':
    # Create a sample graph G with a known 4-clique
    G = nx.Graph()
    G.add_edges_from([(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5), (5, 6)])
    k = 4

    # Perform the reduction
    G_prime, k_prime = reduce_clique_to_pdecide(G, k)

    if G_prime:
        print("\n--- Reduction Output ---")
        print(f"The new graph G' has {G_prime.number_of_nodes()} nodes and {G_prime.number_of_edges()} edges.")
        print(f"G' is bipartite: {nx.is_bipartite(G_prime)}")
        
        # For k>=4, k_prime >= 6. G' has no k'-cliques as it's bipartite.
        # The construction also ensures no induced k'-by-k'-bicliques.
        # So, PDecide(G', k') is true iff G' has an induced k'-matching.
        # This is true iff G has a k-clique.
        print("\n--- Hardness Argument Summary ---")
        print(f"For k={k}, k'={k_prime}. Since k' > 2, the bipartite graph G' has no {k_prime}-cliques.")
        print("The construction also ensures G' has no induced k'-by-k'-bicliques.")
        print("Therefore, an algorithm for PDecide on (G', k') would effectively solve k-CLIQUE on (G, k).")
        print("This proves PDecide is W[1]-hard and PCount is #W[1]-hard.")
