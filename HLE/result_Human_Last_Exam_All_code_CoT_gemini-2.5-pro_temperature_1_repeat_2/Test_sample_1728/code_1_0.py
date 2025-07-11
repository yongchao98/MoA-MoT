def solve_maximal_chromatic_number():
    """
    Calculates and explains the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    """
    # Step 1: Explain the chromatic number property for graph sums.
    # For G = G1 + G2, χ(G) = χ(G1) + χ(G2).
    # For G = C_n + C_n + C_n, this means χ(G) = χ(C_n) + χ(C_n) + χ(C_n).
    print("The chromatic number of a sum of graphs is the sum of their individual chromatic numbers.")
    print("Therefore, the chromatic number of G = C_n + C_n + C_n is χ(G) = χ(C_n) + χ(C_n) + χ(C_n).\n")

    # Step 2: Determine the chromatic number of a cycle C_n.
    # χ(C_n) = 2 for even n, and 3 for odd n.
    print("The chromatic number of a cycle graph, χ(C_n), is:")
    print("- 2, if n is an even number.")
    print("- 3, if n is an odd number.\n")

    # Step 3: Find the maximal chromatic number.
    # To maximize χ(G), we must maximize χ(C_n).
    print("To find the maximal chromatic number of G, we must choose n such that χ(C_n) is maximized.")
    max_chi_cn = 3
    print(f"The maximum value for χ(C_n) is {max_chi_cn}, which occurs when n is odd.\n")

    # Step 4: Calculate the final result and print the equation.
    print("The maximal chromatic number of G is calculated by substituting this value into the equation:")
    result = max_chi_cn + max_chi_cn + max_chi_cn
    print(f"{max_chi_cn} + {max_chi_cn} + {max_chi_cn} = {result}")

solve_maximal_chromatic_number()