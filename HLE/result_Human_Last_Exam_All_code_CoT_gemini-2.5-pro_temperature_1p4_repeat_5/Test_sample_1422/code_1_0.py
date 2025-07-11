def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for a phi^4 interacting scalar field theory.
    """
    # Number of loops given in the problem
    L = 2

    print("To find the minimum number of vertices (V) for a two-loop diagram, we proceed as follows:")
    print("1. Start with the topological formula for loops in a connected diagram: L = I - V + 1")
    print("   (where I = internal lines, V = vertices)")
    print("\n2. For a phi^4 vacuum diagram, the number of internal lines is related to vertices by: I = 2 * V")
    print("\n3. Substituting the second formula into the first gives: L = (2 * V) - V + 1")
    print("   This simplifies to a direct relation between loops and vertices: L = V + 1")
    print("-" * 30)
    
    # We are given L = 2. We solve for V.
    # From L = V + 1, we get V = L - 1
    V = L - 1
    
    print(f"Given the number of loops L = {L}, we solve for V:")
    print(f"Equation: {L} = V + 1")
    print(f"V = {L} - 1")
    print(f"V = {V}")
    print("-" * 30)
    
    # The final equation with all numbers plugged in
    print("The final relationship with the values substituted is:")
    print(f"{L} = {V} + 1")
    
    print("\nTherefore, the minimum number of vertices in a two-loop Feynman diagram for this theory is 1.")
    print("This corresponds to a single-vertex 'figure-eight' vacuum diagram.")

solve_feynman_vertices()