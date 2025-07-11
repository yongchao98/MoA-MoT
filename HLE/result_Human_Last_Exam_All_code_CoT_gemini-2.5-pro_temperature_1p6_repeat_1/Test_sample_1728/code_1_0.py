def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G which is the
    sum of three cycles of length n.
    """

    # The chromatic number of a cycle C_n depends on the parity of n.
    # Case 1: n is even. The chromatic number of an even cycle is 2.
    chi_Cn_even = 2
    # Case 2: n is odd. The chromatic number of an odd cycle is 3.
    chi_Cn_odd = 3

    # The chromatic number of a sum of graphs is the sum of their individual
    # chromatic numbers. So, chi(G) = 3 * chi(C_n).

    # Calculate chi(G) for both cases.
    chi_G_even = 3 * chi_Cn_even
    chi_G_odd = 3 * chi_Cn_odd

    # The maximal chromatic number is the maximum of the two cases.
    maximal_chi_G = max(chi_G_even, chi_G_odd)
    
    # --- Output the explanation ---
    
    print("Let chi(G) be the chromatic number of a graph G.")
    print("The graph G is the sum of three cycles of length n (denoted C_n).")
    print("The chromatic number of a sum of graphs is the sum of their individual chromatic numbers.")
    print("Therefore, chi(G) = chi(C_n) + chi(C_n) + chi(C_n) = 3 * chi(C_n).\n")
    
    print("The value of chi(C_n) depends on whether n is even or odd:")
    
    # Explain the even case
    print(f"1. If n is even, a cycle C_n can be colored with {chi_Cn_even} colors.")
    print(f"   In this case, chi(G) = 3 * {chi_Cn_even} = {chi_G_even}.\n")

    # Explain the odd case
    print(f"2. If n is odd, a cycle C_n requires {chi_Cn_odd} colors.")
    print(f"   In this case, chi(G) = 3 * {chi_Cn_odd} = {chi_G_odd}.\n")

    print(f"To find the maximal chromatic number, we compare the two results: {chi_G_even} and {chi_G_odd}.")
    print(f"The maximum value is {maximal_chi_G}.")
    
    print("\nThe final equation for the maximal case (which occurs when n is odd) is:")
    
    # Final equation showing each number
    num_cycles = 3
    final_chi_cn = chi_Cn_odd
    print(f"chi(G) = {num_cycles} * {final_chi_cn} = {maximal_chi_G}")

solve_chromatic_number()
<<<9>>>