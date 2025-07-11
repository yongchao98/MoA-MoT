def solve_alon_tarsi_knn():
    """
    Calculates the Alon-Tarsi number for the complete bipartite graph K_n,n.
    """
    # Step 1: Identify the value of n for the graph K_{1000,1000}.
    n = 1000

    # Step 2: State the theorem for the Alon-Tarsi number of K_n,n.
    # The theorem is AT(K_n,n) = n + 1.
    print(f"The Alon-Tarsi number of a complete bipartite graph K_n,n is given by the formula:")
    print("AT(K_n,n) = n + 1")
    print("-" * 20)

    # Step 3: Apply the formula for n = 1000.
    alon_tarsi_number = n + 1

    # Step 4: Print the final calculation and the result.
    print(f"For the graph K_{{{n},{n}}}, the value of n is {n}.")
    print("The calculation is as follows:")
    # The final equation with each number is printed below.
    print(f"AT(K_{{{n},{n}}}) = {n} + 1 = {alon_tarsi_number}")

solve_alon_tarsi_knn()