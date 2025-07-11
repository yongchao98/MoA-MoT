def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the
    sum of three cycles of length n.
    """

    # The problem asks for the maximal chromatic number. This occurs when the
    # chromatic number of the cycle, chi(C_n), is maximized.

    # The chromatic number of a cycle C_n is:
    # 2 if n is even (n >= 4)
    # 3 if n is odd (n >= 3)
    # The maximal value is 3.
    max_chi_cn = 3

    # The number of cycles being summed.
    num_cycles = 3

    # The chromatic number of a sum of graphs (G1 + G2 + ...) is the sum of
    # their individual chromatic numbers (chi(G1) + chi(G2) + ...).
    # So, chi(G) = chi(C_n) + chi(C_n) + chi(C_n).
    # The maximal chi(G) is therefore 3 + 3 + 3.
    maximal_chi_g = max_chi_cn * num_cycles

    # Print the explanation and the final equation.
    print(f"The maximal chromatic number of a cycle of length n, chi(C_n), is {max_chi_cn} (when n is odd).")
    print(f"The graph G is the sum of {num_cycles} such cycles.")
    print("The chromatic number of the sum of graphs is the sum of their individual chromatic numbers.")
    print("Therefore, the maximal chromatic number of G is calculated as:")
    print(f"{max_chi_cn} + {max_chi_cn} + {max_chi_cn} = {maximal_chi_g}")

solve_chromatic_number()