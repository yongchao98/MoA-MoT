def solve_maximal_chromatic_number():
    """
    This script calculates the maximal chromatic number of a graph G that is the
    sum (join) of three cycles of length n. It explains the underlying graph
    theory principles and shows the calculation.
    """
    print("Step 1: Understanding the Graph G")
    print("The graph G is the 'sum' of three cycles of length n. In graph theory, this corresponds to the graph join operation.")
    print("So, G = C_n + C_n + C_n, where C_n is a cycle graph of length n.")
    print("-" * 30)

    print("Step 2: Chromatic Number of a Graph Join")
    print("The chromatic number of a join of graphs is the sum of their individual chromatic numbers.")
    print("The formula is: χ(G1 + G2 + ... + Gk) = χ(G1) + χ(G2) + ... + χ(Gk)")
    print("Applying this to our graph G, we get:")
    print("χ(G) = χ(C_n) + χ(C_n) + χ(C_n) = 3 * χ(C_n)")
    print("-" * 30)

    print("Step 3: Chromatic Number of a Cycle C_n")
    print("The chromatic number of a cycle graph C_n depends on whether n is even or odd (assuming n >= 3).")
    
    # Case for even n
    chi_cn_even = 2
    print(f"If n is even, C_n is bipartite, so it can be colored with 2 colors. χ(C_n) = {chi_cn_even}.")
    
    # Case for odd n
    chi_cn_odd = 3
    print(f"If n is odd, C_n requires a third color. χ(C_n) = {chi_cn_odd}.")
    print("-" * 30)

    print("Step 4: Calculating χ(G) for both cases")
    # Calculation for even n
    chi_g_even = 3 * chi_cn_even
    print(f"For even n: χ(G) = 3 * χ(C_n) = 3 * {chi_cn_even} = {chi_g_even}")
    
    # Calculation for odd n
    chi_g_odd = 3 * chi_cn_odd
    print(f"For odd n:  χ(G) = 3 * χ(C_n) = 3 * {chi_cn_odd} = {chi_g_odd}")
    print("-" * 30)
    
    print("Step 5: Finding the Maximal Chromatic Number")
    maximal_chi = max(chi_g_even, chi_g_odd)
    print(f"Comparing the two possible values, {chi_g_even} and {chi_g_odd}, the maximal value is {maximal_chi}.")
    print("This maximum is achieved when n is an odd number.")
    print("\nThe final equation for the maximal case is:")
    print(f"3 * {chi_cn_odd} = {maximal_chi}")

solve_maximal_chromatic_number()