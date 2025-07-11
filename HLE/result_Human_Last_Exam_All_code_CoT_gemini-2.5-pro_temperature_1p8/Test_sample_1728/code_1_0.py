def solve_chromatic_number():
    """
    This function calculates the maximal chromatic number of a graph G,
    where G is the sum of three cycles of length n.

    The problem can be stated as finding the maximal value of chi(C_n + C_n + C_n).

    The chromatic number of a graph join is the sum of the chromatic numbers of the
    individual graphs. Therefore, chi(C_n + C_n + C_n) = chi(C_n) + chi(C_n) + chi(C_n).

    The chromatic number of a cycle, chi(C_n) for n >= 3, is:
    - 2 if n is even (the cycle is bipartite)
    - 3 if n is odd

    To find the maximal chromatic number of G, we need to use the maximal possible
    value for chi(C_n), which is 3. This occurs when n is odd.
    """

    # Maximal chromatic number of a single cycle C_n (when n is odd)
    max_chi_cycle = 3

    # The maximal chromatic number of G is the sum of the maximal chromatic
    # numbers of the three constituent cycles.
    maximal_chi_G = max_chi_cycle + max_chi_cycle + max_chi_cycle

    print("To find the maximal chromatic number of G = C_n + C_n + C_n, we must maximize chi(C_n).")
    print(f"The maximal value for chi(C_n) is {max_chi_cycle}, which occurs when n is odd.")
    print("Therefore, the maximal chromatic number of G is given by the sum:")
    # Print the final equation with each number explicitly shown.
    print(f"{max_chi_cycle} + {max_chi_cycle} + {max_chi_cycle} = {maximal_chi_G}")

solve_chromatic_number()