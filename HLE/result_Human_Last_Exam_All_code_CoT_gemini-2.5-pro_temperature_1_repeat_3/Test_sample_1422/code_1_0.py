def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for an interacting scalar field theory.
    """
    # L = number of loops
    # V = number of vertices
    # I = number of internal lines
    # n = number of lines at each vertex (interaction type)
    
    L = 2
    print(f"The problem states the diagram has L = {L} loops.")
    print("-" * 20)

    # From the topological formula L = I - V + 1, we get I = V + L - 1.
    # For L=2, this is I = V + 1.
    print("Step 1: Use the loop formula L = I - V + 1.")
    print(f"For L={L}, we have 2 = I - V + 1, which rearranges to: I = V + 1.")
    print("-" * 20)

    # For a vacuum diagram (no external lines), n * V = 2 * I.
    # Substituting I = V + 1 gives: n * V = 2 * (V + 1)
    # This simplifies to V * (n - 2) = 2.
    print("Step 2: Use the vertex-line formula n * V = 2 * I (for vacuum diagrams).")
    print("Substituting 'I' from Step 1, we get n * V = 2 * (V + 1).")
    print("Solving for V gives the final equation: V * (n - 2) = 2.")
    print("-" * 20)

    print("Step 3: Find the minimum integer V by testing valid interaction types 'n'.")
    print("An interaction vertex must have at least n=3 lines.")
    
    min_V = float('inf')
    result_details = {}

    # We need V to be a positive integer.
    # From V = 2 / (n - 2), this means (n - 2) must be a positive divisor of 2.
    # The divisors of 2 are 1 and 2.
    for divisor in [1, 2]:
        n_minus_2 = divisor
        n = n_minus_2 + 2
        V = 2 // n_minus_2
        
        print(f"\nTesting case where (n - 2) = {n_minus_2}:")
        print(f"This corresponds to an interaction type n = {n}.")
        print(f"Number of vertices V = 2 / {n_minus_2} = {V}.")
        
        if V < min_V:
            min_V = V
            result_details = {'V': V, 'n': n, 'n_minus_2': n_minus_2}

    print("-" * 20)
    print(f"The minimum number of vertices found is {min_V}.")
    print(f"This is achievable in a phi^{result_details['n']} theory.")
    
    print("\nFinal equation for the minimal case:")
    V_final = result_details['V']
    n_final = result_details['n']
    n_minus_2_final = result_details['n_minus_2']
    
    # Printing each number in the final equation: V * (n - 2) = 2
    print(f"{V_final} * ({n_final} - 2) = 2")
    print(f"{V_final} * {n_minus_2_final} = 2")

solve_feynman_vertices()