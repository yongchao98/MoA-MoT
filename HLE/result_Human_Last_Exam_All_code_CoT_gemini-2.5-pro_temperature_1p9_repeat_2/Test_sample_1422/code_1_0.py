def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for an interacting scalar phi^4 theory.
    """

    # We start with two fundamental formulas for Feynman diagrams:
    # 1. Loop formula: L = I - V + 1
    # 2. Vertex-line relation (for phi^4 theory): 4V = E + 2I
    # where L = loops, I = internal lines, V = vertices, E = external lines.

    print("Step 1: Define the system of equations for a phi^4 theory.")
    print("  L = I - V + 1   (Eq. 1)")
    print("  4V = E + 2I     (Eq. 2)\n")

    # We are given that the number of loops is 2.
    L = 2
    print(f"Step 2: Substitute the given number of loops, L = {L}, into Eq. 1.")
    print(f"  {L} = I - V + 1")
    print("  Solving for I, we get: I = V + 1\n")

    print("Step 3: Substitute the expression for I into Eq. 2.")
    print("  4V = E + 2 * (V + 1)")
    print("  4V = E + 2V + 2\n")

    print("Step 4: Solve the resulting equation for V in terms of E.")
    print("  4V - 2V = E + 2")
    print("  2V = E + 2")
    print("  V = (E + 2) / 2\n")

    # To find the minimum number of vertices (V), we must use the minimum
    # possible number of external lines (E). The minimum is E = 0, which
    # corresponds to a vacuum bubble diagram (a diagram with no external particles).
    E_min = 0
    print(f"Step 5: To find the minimum V, we must use the minimum possible number of external lines, E = {E_min}.")
    
    # Calculate the minimum number of vertices.
    V_min = (E_min + 2) / 2
    
    print("\nFinal Calculation:")
    print(f"The minimum number of vertices is found when E = {E_min}.")
    print(f"V_min = ({E_min} + 2) / 2 = {int(V_min)}")
    print("\nThis corresponds to the 'figure-eight' vacuum diagram, which has one vertex and two loops.")

solve_feynman_vertices()