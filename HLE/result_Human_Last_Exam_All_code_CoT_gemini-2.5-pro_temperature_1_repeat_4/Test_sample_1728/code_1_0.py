def solve_maximal_chromatic_number():
    """
    This function determines the maximal chromatic number of a graph G,
    where G is the sum of three cycles of length n.
    """
    print("The problem is to find the maximal chromatic number of a graph G, where G is the sum of three cycles of length n.")
    print("Let the three cycles be C_n¹, C_n², and C_n³.")
    print("The graph G is their sum: G = C_n¹ + C_n² + C_n³.")
    print("\nStep 1: Use the property of the chromatic number for a sum of graphs.")
    print("The chromatic number of a sum of graphs is the sum of their individual chromatic numbers.")
    print("So, χ(G) = χ(C_n¹) + χ(C_n²) + χ(C_n³) = 3 * χ(C_n).")
    
    print("\nStep 2: Determine the possible values for the chromatic number of a cycle, χ(C_n).")
    print("- If n is even (n >= 4), χ(C_n) = 2.")
    print("- If n is odd (n >= 3), χ(C_n) = 3.")
    
    print("\nStep 3: Find the maximal value of χ(C_n).")
    print("The maximum value for χ(C_n) is 3, which occurs when n is an odd integer.")
    max_chi_cn = 3
    
    print("\nStep 4: Calculate the maximal chromatic number of G.")
    print("To get the maximal χ(G), we use the maximal value of χ(C_n) in our formula.")
    
    num_cycles = 3
    maximal_chromatic_number = num_cycles * max_chi_cn
    
    print("\nThe final calculation is:")
    print(f"{num_cycles} * {max_chi_cn} = {maximal_chromatic_number}")

solve_maximal_chromatic_number()