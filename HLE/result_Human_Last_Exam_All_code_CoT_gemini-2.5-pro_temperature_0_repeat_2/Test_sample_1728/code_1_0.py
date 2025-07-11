def solve_chromatic_number():
    """
    Calculates and explains the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    """
    # Step 1: Explain the relationship between the chromatic number of the sum of graphs
    # and the individual graphs.
    print("The graph G is the sum (join) of three cycles of length n.")
    print("The chromatic number of a graph join is the sum of the chromatic numbers of the individual graphs.")
    print("Therefore, chi(G) = chi(C_n) + chi(C_n) + chi(C_n).")
    print("-" * 30)

    # Step 2: Determine the maximal chromatic number of a single cycle C_n.
    # chi(C_n) = 2 for even n, and 3 for odd n. The maximum is 3.
    max_chi_cycle = 3
    print(f"To find the maximal chromatic number of G, we must maximize chi(C_n).")
    print(f"The chromatic number of a cycle C_n is 2 if n is even, and 3 if n is odd.")
    print(f"The maximum possible chromatic number for a single cycle is {max_chi_cycle}.")
    print("-" * 30)

    # Step 3: Calculate the maximal chromatic number of G.
    num_cycles = 3
    max_chi_G = max_chi_cycle * num_cycles

    # Step 4: Print the final calculation as an equation.
    print("The maximal chromatic number of G is the sum of the maximal chromatic numbers of the three cycles.")
    print(f"Final Equation: {max_chi_cycle} + {max_chi_cycle} + {max_chi_cycle} = {max_chi_G}")

solve_chromatic_number()
<<<9>>>