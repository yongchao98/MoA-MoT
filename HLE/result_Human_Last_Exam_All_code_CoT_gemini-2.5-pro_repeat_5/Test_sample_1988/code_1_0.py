import networkx as nx

def count_induced_subgraphs():
    """
    This script counts the number of induced subgraphs in the Gosset graph
    that are isomorphic to the graph with HoG ID 50698.
    
    The target graph is the complement of the Johnson graph J(8, 2).
    """

    # Step 1: Construct the host graph (Gosset graph) and the target subgraph.
    print("Constructing the Gosset graph (G)...")
    G = nx.gosset_graph()

    print("Constructing the target subgraph (H = complement of Johnson graph J(8,2))...")
    # The vertices of K8 can be represented by the set {0, 1, 2, 3, 4, 5, 6, 7}.
    # The edges of K8 are the 2-element subsets of this set. These are the vertices of H.
    # Two vertices in H are adjacent if the corresponding edges in K8 are disjoint.
    # This is exactly the complement of the Johnson graph J(8, 2).
    H = nx.complement(nx.johnson_graph(8, 2))

    print("-" * 30)
    print(f"Host Graph (G): Gosset graph")
    print(f"  - Vertices: {G.number_of_nodes()}")
    print(f"  - Edges: {G.number_of_edges()}")
    print(f"Target Subgraph (H): Complement of J(8,2)")
    print(f"  - Vertices: {H.number_of_nodes()}")
    print(f"  - Edges: {H.number_of_edges()}")
    print("-" * 30)

    # Step 2: Set up the isomorphism search.
    # We use GraphMatcher, which finds non-induced subgraphs by default.
    # We will add a manual check for the non-edge preservation condition.
    matcher = nx.isomorphism.GraphMatcher(G, H)

    # Pre-calculate the non-edges of H to speed up the check inside the loop.
    non_edges_H = list(nx.non_edges(H))
    
    # This set will store the unique subgraphs we find.
    # A subgraph is identified by its set of vertices.
    # We use frozenset because sets are mutable and cannot be in other sets.
    found_subgraphs = set()

    print("Searching for induced subgraphs. This may take a moment...")
    
    # Step 3: Iterate through all mappings and check for the induced condition.
    for mapping in matcher.subgraph_isomorphisms_iter():
        is_induced = True
        # Check if non-edges in H are preserved as non-edges in G.
        for u, v in non_edges_H:
            # If a non-edge in H maps to an edge in G, it's not an induced subgraph.
            if G.has_edge(mapping[u], mapping[v]):
                is_induced = False
                break
        
        if is_induced:
            # This is a valid induced subgraph isomorphism.
            # Store the set of vertices of this subgraph to count unique ones.
            subgraph_nodes = frozenset(mapping.values())
            found_subgraphs.add(subgraph_nodes)

    # Step 4: Output the final count.
    final_count = len(found_subgraphs)
    
    print("\n" + "=" * 30)
    print("      FINAL RESULT")
    print("=" * 30)
    print("The final count represents the number of unique sets of vertices in the Gosset graph")
    print("that induce a subgraph isomorphic to the target graph.")
    print(f"\nNumber of such subgraphs found: {final_count}")


if __name__ == '__main__':
    count_induced_subgraphs()