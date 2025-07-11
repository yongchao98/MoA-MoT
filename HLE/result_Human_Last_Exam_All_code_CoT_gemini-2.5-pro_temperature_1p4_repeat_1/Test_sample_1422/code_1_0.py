def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in a scalar phi^4 theory.
    """
    print("To find the minimum number of vertices (V) in a two-loop (L=2) Feynman diagram")
    print("for a scalar phi^4 theory, we use two fundamental topological relations:")
    print("1. Loop formula: L = I - V + 1, where I is the number of internal propagators.")
    print("2. Vertex degree formula: 4 * V = 2 * I + E, where E is the number of external lines.")
    print("\nFirst, we combine these equations to express V in terms of L and E.")
    print("From formula 1, we get: I = L + V - 1")
    print("Substituting I into formula 2: 4*V = 2*(L + V - 1) + E")
    print("=> 4*V = 2*L + 2*V - 2 + E")
    print("=> 2*V = 2*L - 2 + E")
    print("=> V = L - 1 + E / 2")
    
    print("\nWe are given L=2. To minimize V, we must choose the smallest possible value for E.")
    print("Since V must be an integer, E must be a non-negative even number.")
    print("The minimum value for E is 0 (a vacuum diagram).")

    # Given values
    L = 2  # Number of loops
    E = 0  # Minimum possible number of external lines

    # Calculate the minimum number of vertices
    V = L - 1 + E / 2
    
    print("\nCalculating the minimum vertices with L=2 and E=0:")
    # The final equation as requested, showing each number
    print(f"Minimum Vertices = {L} - 1 + {E} / 2 = {int(V)}")
    
solve_feynman_vertices()