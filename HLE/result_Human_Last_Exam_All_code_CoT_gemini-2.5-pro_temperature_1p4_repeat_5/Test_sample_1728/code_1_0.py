def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the sum of three cycles of length n.
    """
    
    print("Step 1: Understand the graph G and the 'sum' operation.")
    print("The graph G is the 'sum' of three cycles of length n. This operation is the graph join.")
    print("A key property of the graph join is that the chromatic number of the sum is the sum of the individual chromatic numbers.")
    print("For G = C_n + C_n + C_n, the chromatic number is: χ(G) = χ(C_n) + χ(C_n) + χ(C_n) = 3 * χ(C_n).")
    print("-" * 50)

    print("Step 2: Determine the chromatic number of a cycle C_n.")
    print("The chromatic number of a cycle, χ(C_n), depends on its length n (for n >= 3):")
    print("  - If n is even, χ(C_n) = 2.")
    print("  - If n is odd, χ(C_n) = 3.")
    print("-" * 50)
    
    print("Step 3: Find the maximal chromatic number of G.")
    print("To find the maximal χ(G), we must use the maximal value of χ(C_n).")
    
    # The maximum value for the chromatic number of a cycle C_n is 3.
    chi_cn_max = 3
    print(f"The maximum value for χ(C_n) is {chi_cn_max}, which occurs when n is odd.")
    
    num_cycles = 3
    maximal_chi_g = num_cycles * chi_cn_max
    
    print("\nTherefore, the maximal chromatic number of G is calculated as follows:")
    # Final equation with each number printed, as requested.
    print(f"{num_cycles} * {chi_cn_max} = {maximal_chi_g}")

solve_chromatic_number()