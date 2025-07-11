def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the
    sum (join) of three cycles of length n.
    """

    print("Step 1: Define the graph G and the 'sum' operation.")
    print("The graph G is the sum of three cycles of length n. In graph theory, this 'sum' is commonly interpreted as the graph join operation (+).")
    print("So, G = C_n + C_n + C_n.")
    print("This means G consists of three disjoint n-cycles, and every vertex in each cycle is connected to every vertex in the other two cycles.\n")

    print("Step 2: Use the formula for the chromatic number of a graph join.")
    print("The chromatic number of a join of graphs is the sum of their individual chromatic numbers:")
    print("χ(G) = χ(C_n) + χ(C_n) + χ(C_n)")
    print("χ(G) = 3 * χ(C_n)\n")

    print("Step 3: Determine the chromatic number of a cycle, χ(C_n).")
    print("The chromatic number of a cycle C_n depends on whether n is even or odd:")
    # Chromatic number for an odd cycle (e.g., C_3, C_5, ...)
    chi_cn_odd = 3
    print(f"- If n is odd (n >= 3), χ(C_n) = {chi_cn_odd}.")
    # Chromatic number for an even cycle (e.g., C_4, C_6, ...)
    chi_cn_even = 2
    print(f"- If n is even (n >= 4), χ(C_n) = {chi_cn_even}.\n")

    print("Step 4: Calculate χ(G) for both cases to find the maximum.")
    print("To find the maximal chromatic number of G, we must choose the case for n that maximizes χ(C_n), which is when n is odd.\n")

    # Case 1: n is odd
    num_cycles = 3
    chi_g_odd = num_cycles * chi_cn_odd
    print("Calculation for when n is odd:")
    print(f"χ(G) = {num_cycles} * {chi_cn_odd} = {chi_g_odd}")

    # Case 2: n is even
    chi_g_even = num_cycles * chi_cn_even
    print("\nFor completeness, the calculation for when n is even:")
    print(f"χ(G) = {num_cycles} * {chi_cn_even} = {chi_g_even}\n")

    maximal_chi = max(chi_g_odd, chi_g_even)
    print(f"Comparing the two results ({chi_g_odd} and {chi_g_even}), the maximal chromatic number for G is {maximal_chi}.")


# Run the solver
solve_chromatic_number()