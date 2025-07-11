def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in phi^4 scalar field theory.
    """
    # L: number of loops
    # V: number of vertices
    # I: number of internal lines
    # E: number of external lines

    # The topological relations are:
    # 1) L = I - V + 1
    # 2) 4V = 2I + E (for phi^4 theory)

    # We are given L = 2. We want to find the minimum V.
    L = 2

    # From equation 1, we can express I in terms of L and V:
    # I = L + V - 1
    # Substituting L = 2:
    # I = 2 + V - 1  =>  I = V + 1

    # Now substitute this expression for I into equation 2:
    # 4V = 2 * (V + 1) + E
    # 4V = 2V + 2 + E
    # 2V = 2 + E
    # V = 1 + E / 2

    print("The relationship between vertices (V) and external lines (E) for a two-loop diagram is:")
    print("V = 1 + E / 2")
    print("\nTo find the minimum number of vertices (V), we must find the minimum possible number of external lines (E).")
    print("The number of external lines E must be a non-negative even integer (0, 2, 4, ...).")

    # The minimum possible value for E is 0.
    min_E = 0
    print(f"The minimum valid value for E is {min_E}.")

    # Calculate the minimum V using this minimum E.
    min_V = 1 + min_E / 2

    print("\nSubstituting the minimum E into the equation:")
    print(f"V_min = 1 + {min_E} / 2")
    print(f"V_min = {int(min_V)}")
    print(f"\nThe minimum number of vertices in a two-loop Feynman diagram for this theory is {int(min_V)}.")

solve_feynman_vertices()