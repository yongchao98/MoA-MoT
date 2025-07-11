def solve_alon_tarsi_k1000_1000():
    """
    Calculates the Alon-Tarsi number for the complete bipartite graph K_1000,1000.
    """
    # The graph is K_n,n, where n is 1000.
    n = 1000

    # The Alon-Tarsi number for a complete bipartite graph K_n,n is given by the formula: AT(K_n,n) = n + 1.
    # We will apply this formula to find the answer.
    
    # Calculate the result
    result = n + 1

    # Print the explanation and the final equation.
    print(f"The Alon-Tarsi number for a complete bipartite graph K_n,n is n + 1.")
    print(f"For the graph K_{n},{n} with n = {n}, the calculation is:")
    print(f"{n} + 1 = {result}")

solve_alon_tarsi_k1000_1000()