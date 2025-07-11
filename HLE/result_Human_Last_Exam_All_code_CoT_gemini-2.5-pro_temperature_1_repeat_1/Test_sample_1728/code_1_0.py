def solve_maximal_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the
    sum of three cycles of length n.

    The logic is based on the properties of graph sums and cycle graphs.
    """

    # The chromatic number of a graph sum (join) G1 + G2 is chi(G1) + chi(G2).
    # For G = C_n + C_n + C_n, we have chi(G) = 3 * chi(C_n).

    # The chromatic number of a cycle C_n is:
    # - 2, if n is even (n >= 4)
    # - 3, if n is odd (n >= 3)

    # To find the maximal chromatic number of G, we need the maximal
    # value for the chromatic number of a single cycle, chi(C_n).
    max_chi_cn = 3  # This occurs when n is odd.

    # The number of cycles in the sum.
    num_cycles = 3

    # The maximal chromatic number of G is the sum of the maximal
    # chromatic numbers of the individual cycles.
    maximal_chi_g = num_cycles * max_chi_cn

    print(f"The chromatic number of the sum of three cycles is chi(G) = chi(C_n) + chi(C_n) + chi(C_n).")
    print(f"To find the maximal chromatic number, we need the maximal value for chi(C_n).")
    print(f"The maximal chromatic number for a cycle C_n is {max_chi_cn}, which occurs when n is odd.")
    print(f"\nThus, the maximal chromatic number for G is:")
    print(f"{max_chi_cn} + {max_chi_cn} + {max_chi_cn} = {maximal_chi_g}")


solve_maximal_chromatic_number()
<<<9>>>