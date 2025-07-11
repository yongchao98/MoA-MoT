def solve_maximal_chromatic_number():
    """
    This script explains and calculates the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    """

    # Step 1: Explain the graph structure and the formula for its chromatic number.
    print("Let G be the graph formed by the sum of three cycles of length n, denoted as C_n.")
    print("The 'sum' of graphs corresponds to the graph join operation.")
    print("A key property of the graph join is that the chromatic number of the resulting graph is the sum of the chromatic numbers of the individual graphs.")
    print("Therefore, the chromatic number of G is given by the formula:")
    print("chi(G) = chi(C_n) + chi(C_n) + chi(C_n) = 3 * chi(C_n)\n")

    # Step 2: Analyze the chromatic number of a cycle C_n.
    print("To find the maximal chromatic number of G, we need to find the maximal value of chi(C_n).")
    print("The chromatic number of a cycle C_n (for n >= 3) depends on the parity of n:")
    print(" - If n is even, chi(C_n) = 2 (as the graph is bipartite).")
    print(" - If n is odd, chi(C_n) = 3.\n")

    # Step 3: Determine the maximum value for chi(C_n).
    max_chi_cn = 3
    print(f"The maximum value for the chromatic number of a single cycle C_n is {max_chi_cn}, which occurs when n is an odd number.\n")

    # Step 4: Calculate and display the final answer.
    num_cycles = 3
    max_chi_g = num_cycles * max_chi_cn
    print("Substituting this maximum value into our formula for chi(G), we get the maximal chromatic number:")
    print(f"Maximal chi(G) = {num_cycles} * {max_chi_cn} = {max_chi_g}")


solve_maximal_chromatic_number()
<<<9>>>