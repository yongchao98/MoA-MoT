def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the sum of three cycles of length n.
    """
    print("The problem is to find the maximal chromatic number of a graph G, where G is the sum of three cycles of length n.")
    print("Let the three cycles be C_n^1, C_n^2, and C_n^3. Then G = C_n^1 + C_n^2 + C_n^3.")
    print("\nStep 1: Understand the chromatic number of a graph sum (join).")
    print("The chromatic number of a join of graphs is the sum of their individual chromatic numbers.")
    print("So, X(G) = X(C_n^1) + X(C_n^2) + X(C_n^3).")
    
    num_cycles = 3
    print(f"Since all three graphs are cycles of length n, X(G) = {num_cycles} * X(C_n).\n")
    
    print("Step 2: Understand the chromatic number of a cycle C_n (for n >= 3).")
    print("- If n is even, X(C_n) = 2.")
    print("- If n is odd, X(C_n) = 3.")
    
    chi_cn_even = 2
    chi_cn_odd = 3
    
    print("\nStep 3: Find the maximal chromatic number of G.")
    print("To maximize X(G), we must maximize X(C_n).")
    max_chi_cn = max(chi_cn_even, chi_cn_odd)
    print(f"The maximum value for X(C_n) is {max_chi_cn}, which occurs when n is odd.")
    
    print("\nStep 4: Calculate the final result.")
    maximal_chromatic_number = num_cycles * max_chi_cn
    
    # Final output showing the numbers in the equation
    print("The maximal chromatic number of G is calculated as:")
    print(f"{num_cycles} * {max_chi_cn} = {maximal_chromatic_number}")

solve_chromatic_number()
<<<9>>>