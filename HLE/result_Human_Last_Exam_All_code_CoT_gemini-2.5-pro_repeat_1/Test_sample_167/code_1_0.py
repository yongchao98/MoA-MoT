def solve_alon_tarsi_k_mn():
    """
    Calculates the Alon-Tarsi number for the complete bipartite graph K_1000,1000.
    """
    # The graph is K_m,n
    m = 1000
    n = 1000

    # For a complete bipartite graph K_m,n, the Alon-Tarsi number is min(m, n) + 1.
    min_partition_size = min(m, n)
    alon_tarsi_number = min_partition_size + 1

    # Print the equation and the final result.
    print(f"The Alon-Tarsi number for K_{m},{n} is given by the formula: min(m, n) + 1")
    print(f"min({m}, {n}) + 1 = {alon_tarsi_number}")

solve_alon_tarsi_k_mn()