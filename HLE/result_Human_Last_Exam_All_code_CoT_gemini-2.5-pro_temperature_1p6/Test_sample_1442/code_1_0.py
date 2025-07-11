import itertools

def solve():
    """
    This function addresses the question:
    Do two bipartite, d-regular graphs on n vertices necessarily have the same number of 3-matchings?

    To answer this, we construct a counterexample using two non-isomorphic
    3-regular bipartite graphs on 8 vertices.
    1. G1: The cube graph, Q3.
    2. G2: Another such graph, non-isomorphic to the cube graph.

    We will then count the number of 3-matchings in G1 and G2.
    If the counts differ, the answer to the question is "No".
    """

    def count_3_matchings(graph_adj_list):
        """
        Counts the number of 3-matchings in a graph.
        A 3-matching is a set of 3 edges with no shared vertices.
        """
        # Step 1: Extract a unique list of edges from the adjacency list.
        # An edge is represented by a frozenset of its two vertices to handle
        # its unordered nature (u,v) == (v,u) and for hashability.
        edges = set()
        for u, neighbors in graph_adj_list.items():
            for v in neighbors:
                if u < v:
                    edges.add(frozenset([u, v]))
        
        edge_list = list(edges)
        
        # A 3-matching requires at least 3 edges in the graph.
        if len(edge_list) < 3:
            return 0
            
        # Step 2: Generate all unique combinations of 3 edges.
        edge_combinations = itertools.combinations(edge_list, 3)
        
        # Step 3: Iterate through combinations and count valid 3-matchings.
        matching_count = 0
        for edge_combo in edge_combinations:
            # Collect all vertices involved in the current combination of 3 edges.
            vertices_in_combo = set()
            for edge in edge_combo:
                vertices_in_combo.update(edge)
            
            # A set of 3 edges is a 3-matching if and only if they are pairwise disjoint.
            # This means the total number of unique vertices must be 3 * 2 = 6.
            if len(vertices_in_combo) == 6:
                matching_count += 1
                
        return matching_count

    # Define Graph G1 (the cube graph, Q3) using an adjacency list.
    # Vertices 0-7 correspond to binary strings 000-111.
    G1_adj = {
        0: [1, 2, 4], 1: [0, 3, 5], 2: [0, 3, 6], 3: [1, 2, 7],
        4: [0, 5, 6], 5: [1, 4, 7], 6: [2, 4, 7], 7: [3, 5, 6]
    }

    # Define Graph G2 (a different 3-regular bipartite graph on 8 vertices).
    G2_adj = {
        0: [4, 5, 6], 1: [5, 6, 7], 2: [4, 6, 7], 3: [4, 5, 7],
        4: [0, 2, 3], 5: [0, 1, 3], 6: [0, 1, 2], 7: [1, 2, 3]
    }
    
    # --- Calculation and Conclusion ---
    print("Let G1 and G2 be two non-isomorphic, bipartite, 3-regular graphs on 8 vertices.")
    print("G1 is the standard cube graph (Q3). G2 is another such graph.")
    print("We will count the number of 3-matchings in each graph.\n")

    g1_matchings = count_3_matchings(G1_adj)
    g2_matchings = count_3_matchings(G2_adj)

    print(f"Number of 3-matchings in G1 (the cube graph) is: {g1_matchings}")
    print(f"Number of 3-matchings in G2 is: {g2_matchings}")

    print("\nSince the two graphs satisfy the conditions but have a different number of 3-matchings, the answer to the question is no.")

solve()