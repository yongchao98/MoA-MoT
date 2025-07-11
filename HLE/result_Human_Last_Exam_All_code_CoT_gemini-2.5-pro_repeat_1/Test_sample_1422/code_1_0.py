def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in an interacting scalar field theory (phi^4).
    """

    # Number of loops required
    L = 2

    # To find the minimum number of vertices (V), we must use the minimum
    # possible number of external lines (E). The absolute minimum is E=0,
    # which corresponds to a vacuum diagram.
    E = 0

    # The general formula relating Vertices (V), Loops (L), and
    # External lines (E) for a phi^4 theory is: V = L - 1 + E / 2
    # We will now calculate V using the values for L and E.

    min_V = L - 1 + E / 2

    print("Problem: Find the minimum number of vertices in a two-loop Feynman diagram for a phi^4 scalar theory.")
    print("--------------------------------------------------------------------------------------------")
    print("The relationship between Vertices (V), Loops (L), and External lines (E) is given by the formula:")
    print("V = L - 1 + E / 2")
    print("\nWe are given L=2. To minimize V, we must minimize E.")
    print("The minimum possible number of external lines is E=0.")
    print("\nSubstituting L=2 and E=0 into the formula:")
    # We use int() because variables in the equation are integers
    print(f"V = {L} - 1 + {int(E)} / 2")
    print(f"V = {int(min_V)}")
    print("\nTherefore, the minimum number of vertices is 1.")
    print("This corresponds to a 'figure-eight' vacuum diagram where a single vertex's four lines connect back to itself.")

solve_feynman_vertices()