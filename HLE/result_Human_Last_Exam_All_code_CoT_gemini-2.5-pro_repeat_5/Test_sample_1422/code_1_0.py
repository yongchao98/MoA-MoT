def find_min_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in an interacting scalar field theory (assumed to be phi^4).
    """

    # --- Problem Setup ---
    # L: Number of loops
    # k: Number of lines per vertex (4 for phi^4 theory)
    L_target = 2
    k = 4

    print("Goal: Find the minimum number of vertices (V) for a Feynman diagram with L=2 loops.")
    print("Assumption: The theory is a phi^4 scalar theory, so each vertex connects k=4 lines.")
    print("\nThe topological relations for a connected diagram are:")
    print("1. Loop formula: L = I - V + 1")
    print("2. Vertex-line counting: k*V = 2*I + E")
    print("where I is the number of internal lines and E is the number of external lines.")
    print("-" * 50)
    print("We will search for the smallest integer V >= 1 that satisfies these equations for L=2.")

    # --- Search Algorithm ---
    # Iterate through possible numbers of vertices, starting from the minimum possible, V=1.
    for V in range(1, 10):
        # For each V, iterate through possible numbers of external lines, starting from E=0.
        # A diagram with no external lines (E=0) is a "vacuum diagram".
        for E in range(0, 10):
            # From the vertex-line counting formula: 2*I = k*V - E
            # For I to be a non-negative integer, (k*V - E) must be a non-negative and even number.
            if (k * V - E) >= 0 and (k * V - E) % 2 == 0:
                I = (k * V - E) // 2

                # Now check if these values satisfy the loop formula
                if (I - V + 1) == L_target:
                    print(f"Found a valid configuration with V = {V} vertices.")
                    print(f"This configuration has E = {E} external lines and I = {I} internal lines.")
                    print("\nVerifying the equations:")
                    print("Equation 1 (Loops): L = I - V + 1")
                    print(f"Plugging in the numbers: {L_target} = {I} - {V} + 1")
                    print(f"Result: {L_target} = {I - V + 1}, which is correct.")
                    
                    print("\nEquation 2 (Vertices): k*V = 2*I + E")
                    print(f"Plugging in the numbers: {k}*{V} = 2*{I} + {E}")
                    print(f"Result: {k*V} = {2*I + E}, which is correct.")

                    print(f"\nSince we searched starting from V=1, the minimum number of vertices is {V}.")
                    
                    # Return the final answer
                    return V
    return None

min_v = find_min_vertices()
# The problem asks for the answer in a specific format at the end.
# The code above explains the derivation and prints the result.
# The final answer is the value returned by the function.
print(f"\n<<<1>>>")
