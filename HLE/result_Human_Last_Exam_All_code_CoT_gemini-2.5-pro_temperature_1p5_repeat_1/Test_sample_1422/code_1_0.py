def find_minimum_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    based on the topological formula for Feynman diagrams.
    """
    
    # The number of loops required.
    L = 2
    
    # The equation relating Vertices (V), Loops (L), interaction type (n),
    # and External lines (E) is: (n - 2) * V = 2*L - 2 + E
    
    print("The general formula for a Feynman diagram is: (n - 2) * V = 2*L - 2 + E")
    print(f"We are looking for the minimum number of vertices (V) for a two-loop (L={L}) diagram.\n")
    
    # To minimize V, we must choose appropriate values for n and E.
    # - We minimize the number of external lines, E. The minimum is E=0 (a vacuum diagram).
    # - We test different interaction types, n (where n >= 3).
    
    # Case: Quartic interaction (n=4), like in phi-4 theory, and a vacuum diagram (E=0).
    n = 4
    E = 0
    
    print(f"Let's analyze the case with the minimum external lines (E = {E}) for a theory with a quartic interaction (n = {n}).")

    # Calculate the left and right sides of the simplified equation: (n - 2) * V = 2 + E
    lhs_factor = n - 2
    rhs = 2 + E

    print(f"Plugging in the values, the equation becomes: ({n} - {2}) * V = {2} + {E}")
    print(f"This simplifies to: {lhs_factor} * V = {rhs}")
    
    # Solve for V
    if rhs % lhs_factor == 0:
        V = rhs // lhs_factor
        print(f"Solving for V gives V = {rhs} / {lhs_factor}, which results in V = {V}.")
        print(f"\nThe minimum number of vertices is {V}.")
    else:
        # This part will not be reached for our chosen n and E
        print("This combination of n and E does not result in an integer number of vertices.")

find_minimum_vertices()
<<<1>>>