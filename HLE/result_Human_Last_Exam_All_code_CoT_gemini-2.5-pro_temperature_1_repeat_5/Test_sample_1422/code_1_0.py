def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for an interacting scalar phi^4 field theory.
    """
    # L = number of loops
    # V = number of vertices
    # I = number of internal lines
    # E = number of external lines

    print("To find the minimum number of vertices (V), we use two key formulas for connected Feynman diagrams:")
    print("1. Loop Formula: L = I - V + 1")
    print("2. Vertex Rule (phi^4 theory): 4 * V = 2 * I + E")
    print("-" * 60)

    print("We are given L = 2 (for a two-loop diagram).")
    print("To minimize V, we consider the simplest case with a minimum number of external lines, E = 0 (a vacuum bubble diagram).")
    print("\nOur system of equations is:")
    print("1. 2 = I - V + 1")
    print("2. 4 * V = 2 * I + 0")
    print("-" * 60)

    print("Step 1: Rearrange the first equation to solve for I.")
    print("From '2 = I - V + 1', we get 'I = V + 1'.")
    print("\nStep 2: Substitute 'I = V + 1' into the second equation '4 * V = 2 * I'.")
    # This represents the equation 4*V = 2*(V+1)
    print("This gives us the equation: 4 * V = 2 * (V + 1)")
    print("\nStep 3: Solve the equation for V.")
    # This represents 4*V = 2*V + 2
    print("Expanding the right side: 4 * V = 2 * V + 2")
    # This represents 4*V - 2*V = 2
    print("Subtracting '2 * V' from both sides: 2 * V = 2")
    # This represents V = 1
    v = 1
    print("Dividing by 2: V = {}".format(v))
    print("-" * 60)
    print("This result corresponds to the 'figure-eight' vacuum diagram, which has 1 vertex and 2 loops.")
    print("\nTherefore, the minimum number of vertices is 1.")

solve_feynman_vertices()