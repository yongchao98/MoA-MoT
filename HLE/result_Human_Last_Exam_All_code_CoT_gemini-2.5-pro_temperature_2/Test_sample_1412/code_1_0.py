def solve():
    """
    This problem is about finding the number of non-isomorphic graphs with a specific set of properties.
    The properties are: connected, 3-regular, 2000 vertices, having a perfect matching, and being "adjustable".
    
    1.  A key insight is the characterization of such graphs. A graph with these properties must be a "decorated prism graph".
    2.  The structure of a decorated prism graph on 2n=2000 vertices is determined by a partition of the edges of a cycle C_1000 into two subgraphs, H and G.
    3.  The "adjustable" property imposes very strong local constraints on the edge partition of C_1000. For any edge, its two neighbors in the cycle must have the same color (i.e., both belong to H or both belong to G).
    4.  This local rule forces the global coloring pattern of C_1000 to be one of two types:
        a) Monochromatic: All edges have the same color (all H or all G).
        b) Alternating: The colors alternate around the cycle (H, G, H, G, ...).
    5.  These three distinct coloring schemes lead to three non-isomorphic graphs:
        - The 'all-H' coloring gives the prism graph C_1000 x K_2.
        - The 'all-G' coloring gives a bipartite graph which is non-isomorphic to the prism graph.
        - The alternating coloring gives a third graph, non-isomorphic to the first two.
        
    Therefore, there are exactly 3 such non-isomorphic graphs.
    """
    
    # The number of non-isomorphic graphs resulting from the analysis
    num_graphs = 3
    
    # Each possibility corresponds to a specific graph structure
    # Let's print the number and the reasoning for it.
    
    print("The number of non-isomorphic graphs is based on the possible structures of the 'decorated prism graph'.")
    print("Let the base graph be C_1000. Its edges are partitioned into two sets, H and G.")
    
    # Case 1: All edges are in H (monochromatic)
    graph_1_edges_H = 1000
    graph_1_edges_G = 0
    print(f"1. Monochromatic H coloring: H has {graph_1_edges_H} edges (a C_1000), G has {graph_1_edges_G} edges. This gives the prism graph C_1000 x K_2.")
    
    # Case 2: All edges are in G (monochromatic)
    graph_2_edges_H = 0
    graph_2_edges_G = 1000
    print(f"2. Monochromatic G coloring: H has {graph_2_edges_H} edges, G has {graph_2_edges_G} edges (a C_1000). This gives a non-isomorphic bipartite graph.")
    
    # Case 3: Alternating colors
    graph_3_edges_H = 500
    graph_3_edges_G = 500
    print(f"3. Alternating H/G coloring: H has {graph_3_edges_H} edges (a perfect matching), G has {graph_3_edges_G} edges (a perfect matching). This gives a third, non-isomorphic graph.")
    
    # The equation represents the sum of distinct cases found.
    print(f"Total number of non-isomorphic graphs = 1 (all H) + 1 (all G) + 1 (alternating) = {num_graphs}")


solve()