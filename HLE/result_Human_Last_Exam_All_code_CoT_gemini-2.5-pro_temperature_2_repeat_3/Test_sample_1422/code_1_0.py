def find_minimum_vertices():
    """
    This function searches for the minimum number of vertices (V) required
    to form a two-loop (L=2) Feynman diagram in an interacting scalar field theory.
    It iterates through possible values for V, interaction type (n), and
    external lines (E) to find the first combination that satisfies the
    topological rules of Feynman diagrams.
    """
    
    loops = 2

    # We want to find the minimum V, so we start the search with V=1.
    for v in range(1, 10):  # Test for V = 1, 2, 3...
        
        # An interaction requires at least 3 lines meeting at a vertex.
        for n in range(3, 10): # Test interaction type n=3 (phi^3), n=4 (phi^4)...
            
            # Diagrams can have 0 or more external lines (E=0 is a vacuum bubble).
            for e in range(0, 10): # Test E = 0, 1, 2...
                
                # From the vertex rule: n * V = 2 * I + E
                # Therefore, 2 * I = n * V - E
                # The term (n * V - E) must be a non-negative, even integer
                # for I (number of internal lines) to be a valid non-negative integer.
                
                term = n * v - e
                if term >= 0 and term % 2 == 0:
                    i = term // 2
                    
                    # Now check the topological loop rule: L = I - V + 1
                    if loops == i - v + 1:
                        # We have found a valid set of parameters. Since we
                        # are iterating V from 1 upwards, this first solution
                        # corresponds to the minimum number of vertices.
                        print("Found a valid diagram with the minimum number of vertices.")
                        print("-" * 55)
                        print(f"Number of Vertices (V): {v}")
                        print(f"Number of Loops (L): {loops}")
                        print(f"Interaction Type (n-point): {n}")
                        print(f"Number of External Lines (E): {e}")
                        print(f"Number of Internal Lines (I): {i}")
                        print("\nVerification using the topological formulas:")
                        
                        # Print the final equations with the numbers filled in
                        print("\n1. Vertex Rule: n * V = 2 * I + E")
                        print(f"                {n} * {v} = 2 * {i} + {e}")
                        print(f"                {n*v} = {2*i + e}  (Rule Satisfied)")
                        
                        print("\n2. Loop Rule: L = I - V + 1")
                        print(f"             {loops} = {i} - {v} + 1")
                        print(f"             {loops} = {i - v + 1}  (Rule Satisfied)")
                        
                        # Return the minimum value to be used in the final answer
                        return v
    return None

min_vertices = find_minimum_vertices()

# The final answer is wrapped in <<<>>>
print(f"\n<<<The minimum number of vertices is {min_vertices}.>>>")