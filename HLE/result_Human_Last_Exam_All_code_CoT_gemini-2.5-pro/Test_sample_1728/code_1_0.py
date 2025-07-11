def solve_chromatic_number_sum_of_cycles():
    """
    Calculates the maximal chromatic number of a graph G that is the
    sum of three cycles of length n.
    """

    # The number of cycles being summed.
    num_cycles = 3

    # The chromatic number of a cycle C_n depends on its length n.
    # Case 1: n is even (n >= 4). The cycle is bipartite.
    chi_Cn_even = 2
    # Case 2: n is odd (n >= 3).
    chi_Cn_odd = 3

    # To find the maximal chromatic number of the resulting graph G, we need to
    # find the maximal chromatic number of a single cycle C_n.
    max_chi_Cn = max(chi_Cn_even, chi_Cn_odd)

    # The chromatic number of a sum of graphs (graph join) is the sum of
    # their individual chromatic numbers.
    # χ(G) = χ(C_n) + χ(C_n) + χ(C_n) = num_cycles * χ(C_n)
    # The maximal chromatic number of G is therefore num_cycles * max_chi_Cn.
    max_chi_G = num_cycles * max_chi_Cn

    # Print the explanation and the final equation.
    print(f"The chromatic number of a cycle of length n, χ(C_n), can be {chi_Cn_even} or {chi_Cn_odd}.")
    print(f"The maximal chromatic number for a single cycle, max(χ(C_n)), is {max_chi_Cn}.")
    print("The chromatic number of the sum of three cycles is 3 * χ(C_n).")
    print("Therefore, the maximal chromatic number of G is:")
    print(f"{num_cycles} * {max_chi_Cn} = {max_chi_G}")

solve_chromatic_number_sum_of_cycles()