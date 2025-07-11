def solve_alon_tarsi_kmn():
    """
    Calculates the Alon-Tarsi number for a complete bipartite graph K_m,n.
    """
    # Define the parameters for the complete bipartite graph K_m,n
    m = 1000
    n = 1000

    # The Alon-Tarsi number of K_m,n is given by the formula: min(m, n) + 1.
    min_val = min(m, n)
    result = min_val + 1

    # Print the explanation and the final equation with all numbers.
    print(f"The Alon-Tarsi number of the complete bipartite graph K_{m},{n} is found using the formula min(m, n) + 1.")
    print("For this case, m = 1000 and n = 1000.")
    print("The calculation is:")
    print(f"min({m}, {n}) + 1 = {min_val} + 1 = {result}")

solve_alon_tarsi_kmn()