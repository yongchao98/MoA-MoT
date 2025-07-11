def solve_maximal_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the sum 
    of three cycles of length n.
    """
    
    print("Step 1: Define the graph G and the chromatic number of a graph join.")
    print("The graph G is the sum (join) of three cycles of length n: G = C_n + C_n + C_n.")
    print("The chromatic number of a graph join is the sum of the individual chromatic numbers.")
    print("Therefore, chi(G) = chi(C_n) + chi(C_n) + chi(C_n) = 3 * chi(C_n).")
    print("-" * 30)

    print("Step 2: Determine the possible chromatic numbers for a single cycle C_n.")
    print("The chromatic number of a cycle C_n depends on the parity of n (for n >= 3):")
    print("- chi(C_n) = 2, if n is even.")
    print("- chi(C_n) = 3, if n is odd.")
    print("-" * 30)

    print("Step 3: Find the maximal chromatic number for G.")
    print("To find the maximal value for chi(G), we must use the maximal value for chi(C_n).")
    
    num_cycles = 3
    max_chi_cn = 3 # The maximum chromatic number for a single cycle occurs when n is odd.
    
    print(f"The number of cycles is {num_cycles}.")
    print(f"The maximal chromatic number of a single cycle is {max_chi_cn}.")
    
    maximal_chromatic_number_G = num_cycles * max_chi_cn
    
    print("\nFinal Calculation:")
    print(f"The maximal chromatic number of G is found by the equation:")
    # The final equation with each number printed, as requested.
    print(f"{num_cycles} * {max_chi_cn} = {maximal_chromatic_number_G}")

solve_maximal_chromatic_number()