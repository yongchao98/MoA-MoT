def solve_maximal_chromatic_number():
    """
    This script calculates the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    """

    # The chromatic number of a sum (join) of graphs is the sum of their individual chromatic numbers.
    # For a graph G that is the sum of three cycles of length n (G = C_n + C_n + C_n),
    # the chromatic number is chi(G) = chi(C_n) + chi(C_n) + chi(C_n).

    # To find the maximal chromatic number of G, we need the maximal chromatic number of a cycle C_n.
    # A cycle C_n has a chromatic number of 2 if n is even, and 3 if n is odd (for n >= 3).
    # The maximum of these two values is 3.
    max_chi_cn = 3
    
    # The number of cycles being summed is 3.
    num_cycles = 3
    
    # The maximal chromatic number of G is the sum of the maximal chromatic numbers of the three cycles.
    total_max_chi = max_chi_cn * num_cycles
    
    # Construct the equation string as requested.
    equation_parts = [str(max_chi_cn)] * num_cycles
    equation_str = " + ".join(equation_parts)
    
    print(f"The maximal chromatic number for one cycle C_n is {max_chi_cn}.")
    print("The maximal chromatic number for the sum of three cycles is the sum of their individual maximal chromatic numbers.")
    print("The final calculation is:")
    print(f"{equation_str} = {total_max_chi}")

solve_maximal_chromatic_number()