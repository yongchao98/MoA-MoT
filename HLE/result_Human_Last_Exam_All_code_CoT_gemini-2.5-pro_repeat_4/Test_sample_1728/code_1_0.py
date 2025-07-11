def solve_maximal_chromatic_number():
    """
    Calculates and explains the maximal chromatic number of a graph G,
    where G is the sum of three cycles of length n.
    """
    print("This program calculates the maximal chromatic number of a graph G = C_n + C_n + C_n.")
    print("The '+' operation here denotes the disjunctive product of graphs.")
    print("-" * 60)

    # Step 1: Explain the formula for the chromatic number of a graph sum.
    print("Step 1: Use the property of the chromatic number for a disjunctive product.")
    print("The chromatic number of a sum of graphs is the product of their individual chromatic numbers.")
    print("χ(G1 + G2) = χ(G1) * χ(G2)")
    print("Therefore, for G = C_n + C_n + C_n, we have:")
    print("χ(G) = χ(C_n) * χ(C_n) * χ(C_n) = (χ(C_n))^3")
    print("-" * 60)

    # Step 2: Explain the chromatic number of a cycle C_n.
    print("Step 2: Determine the chromatic number of a cycle C_n.")
    print("The chromatic number of a cycle C_n (for n >= 3) depends on the parity of n:")
    print(" - If n is even, C_n is bipartite, and its chromatic number is 2.")
    print(" - If n is odd, C_n requires three colors, and its chromatic number is 3.")
    print("-" * 60)
    
    # Step 3: Find the maximum value for χ(C_n).
    print("Step 3: Find the maximal value for χ(C_n).")
    print("To maximize χ(G), we must maximize χ(C_n). The maximum value for χ(C_n) is 3.")
    max_chi_cn = 3
    print(f"Maximal χ(C_n) = {max_chi_cn} (when n is odd).")
    print("-" * 60)

    # Step 4: Calculate and display the final result.
    print("Step 4: Calculate the maximal chromatic number for G.")
    num_cycles = 3
    result = max_chi_cn ** num_cycles

    print("The maximal chromatic number of G is (Maximal χ(C_n))^3.")
    # Printing the final equation with each number, as requested.
    equation_parts = [str(max_chi_cn)] * num_cycles
    print(f"Final Calculation: {max_chi_cn} * {max_chi_cn} * {max_chi_cn} = {result}")

# Execute the function
solve_maximal_chromatic_number()