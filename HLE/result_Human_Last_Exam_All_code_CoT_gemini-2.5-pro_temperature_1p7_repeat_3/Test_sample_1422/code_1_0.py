def find_minimum_vertices():
    """
    Finds the minimum number of vertices (V) for a two-loop (L=2) Feynman diagram
    in an interacting scalar field theory (phi^n).

    The search is based on the derived topological relation:
    V * (n - 2) = E + 2
    where:
    V = number of vertices (must be integer >= 1)
    n = order of interaction (must be integer >= 3)
    E = number of external lines (must be integer >= 0)
    """
    # We are looking for the minimum possible integer V, so we start our search at V=1.
    for V_candidate in range(1, 10):
        # For a given V, we check if there exists any valid theory (n>=3) and
        # any valid number of external lines (E>=0) that satisfy the equation.
        # Let's check a few simple interaction types and external line counts.
        for n in range(3, 10):
            for E in range(0, 10):
                # Check if the current combination of (V, n, E) is a valid solution
                if V_candidate * (n - 2) == E + 2:
                    print(f"Solution found: A diagram with V = {V_candidate} vertices is possible.")
                    print(f"This occurs in a phi^{n} theory with E = {E} external lines.")
                    
                    # The number of internal lines is determined by I = V + 1 (since L=2)
                    I = V_candidate + 1
                    L = 2
                    print(f"This {L}-loop diagram would have I = {I} internal lines.")

                    print("\nVerifying with the core equation: V * (n - 2) = E + 2")
                    print(f"Substituting the values: {V_candidate} * ({n} - 2) = {E} + 2")
                    left_side = V_candidate * (n-2)
                    right_side = E + 2
                    print(f"Final check: {left_side} = {right_side}")
                    
                    # Because we are iterating V_candidate from 1 upwards, the first
                    # solution we find gives the minimum possible number of vertices.
                    print(f"\nThe minimum number of vertices is {V_candidate}.")
                    return V_candidate

# Run the function to find and print the result.
min_vertices = find_minimum_vertices()

# Final answer in the required format.
print(f"<<<{min_vertices}>>>")
