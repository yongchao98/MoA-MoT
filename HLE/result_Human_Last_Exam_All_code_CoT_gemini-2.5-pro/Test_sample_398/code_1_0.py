import networkx as nx

def solve_graph_problem():
    """
    Analyzes the properties of a graph G' = G + e, where G is maximal planar.
    """
    # 1. & 2. Create a maximal planar graph (Octahedron) and verify its properties.
    G = nx.octahedral_graph()
    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    print("Step 1 & 2: Analyze a maximal planar graph G (Octahedron)")
    print(f"Number of vertices (n): {n}")
    print(f"Number of edges (m): {m}")
    
    # For a maximal planar graph with n>=3, m = 3n - 6.
    expected_m = 3 * n - 6
    print(f"Checking the condition for a maximal planar graph: m = 3*n - 6")
    print(f"Equation: {m} = 3 * {n} - 6")
    print(f"Result: {m} = {expected_m}")

    is_planar = nx.is_planar(G)
    is_maximal = (m == expected_m)

    if is_planar and is_maximal:
        print("G is confirmed to be a maximal planar graph.")
    else:
        print("The chosen graph is not maximal planar. Aborting.")
        return
    print("-" * 30)

    # 3. Add a non-edge 'e' to create G' and check its planarity.
    # In the nx.octahedral_graph, the non-edges are pairs of opposite vertices.
    # Let's pick e = (0, 5), which connects two opposite apexes.
    # Note: A different choice, like (1, 3), would also work.
    e_to_add = (0, 5) 
    
    # Check if it's a non-edge before adding
    if G.has_edge(*e_to_add):
        # Find the actual non-edges in the complement graph
        non_edges = list(nx.complement(G).edges())
        e_to_add = non_edges[0]

    print(f"Step 3: Add a non-edge e = {e_to_add} to G to form G'")
    G_prime = G.copy()
    G_prime.add_edge(*e_to_add)

    is_g_prime_planar = nx.is_planar(G_prime)
    print(f"Is the new graph G' planar? {is_g_prime_planar}")
    if not is_g_prime_planar:
        print("Conclusion: G' is not planar, so it must be drawn with at least one crossing.")
    else:
        print("Error: G' should be non-planar.")
        return
    print("-" * 30)
    
    # 4. & 5. Check if a 1-crossing drawing exists and if it's unique.
    # This is done by finding which edges e' from G can be "swapped" with e.
    # If G - e' + e is planar, it means G' can be drawn with e and e' crossing.
    possible_crossed_edges = []
    for e_prime in G.edges():
        H = G.copy()
        H.remove_edge(*e_prime)
        H.add_edge(*e_to_add)
        if nx.is_planar(H):
            possible_crossed_edges.append(e_prime)
    
    print("Step 4 & 5: Investigate 1-crossing drawings of G'")
    if len(possible_crossed_edges) > 0:
        print(f"A drawing with one crossing is possible.")
        print(f"The new edge {e_to_add} can be drawn by crossing any of the following {len(possible_crossed_edges)} edges:")
        for edge in possible_crossed_edges:
            print(f"  - {edge}")
        
        if len(possible_crossed_edges) > 1:
            print("\nConclusion: Since there is more than one choice for the crossed edge, the drawing with one crossing is not unique.")
        else:
            print("\nConclusion: The choice of edge to cross is unique.")
    else:
        print("Conclusion: A drawing with a single crossing was not found. This may indicate G' requires more than one crossing.")

    print("-" * 30)
    print("Final Summary:")
    print("1. G' is non-planar (requires at least one crossing).")
    print("2. G' can be drawn with exactly one crossing.")
    print("3. The drawing is not unique because there are multiple choices for which edge to cross.")
    print("This matches statement B.")

solve_graph_problem()