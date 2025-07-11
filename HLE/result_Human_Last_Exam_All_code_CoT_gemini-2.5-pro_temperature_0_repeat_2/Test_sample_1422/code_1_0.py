def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for an interacting scalar field theory.
    """
    # Step 1 & 2: State and combine the topological formulas.
    # L = I - V + 1  =>  I = L + V - 1
    # n * V = 2 * I + E
    # Substituting I into the second equation:
    # n * V = 2 * (L + V - 1) + E
    # n * V = 2*L + 2*V - 2 + E
    # V * (n - 2) = 2*L - 2 + E
    # V = (2*L - 2 + E) / (n - 2)
    print("The general formula relating vertices (V), loops (L), external lines (E),")
    print("and interaction type (n, for a phi^n theory) is:")
    print("V = (2*L - 2 + E) / (n - 2)\n")

    # Step 3: Apply the constraint L=2.
    L = 2
    print(f"For a two-loop diagram, we set L = {L}:")
    print("V = (2*2 - 2 + E) / (n - 2)")
    print("V = (2 + E) / (n - 2)\n")

    # Step 4: Find the minimum integer V.
    # To minimize V, we must minimize the numerator (2 + E) and maximize the
    # denominator (n - 2).
    # The smallest possible value for E (number of external lines) is 0.
    # The interaction type n must be an integer >= 3.
    # The number of vertices V must be a positive integer.
    print("To find the minimum V, we search for the smallest positive integer V")
    print("that can be formed with valid parameters n (>=3) and E (>=0).\n")

    # We can test for V=1, the lowest possible positive integer.
    min_V = 1
    # If V = 1, the equation becomes:
    # 1 = (2 + E) / (n - 2)
    # n - 2 = 2 + E
    # n = 4 + E

    # We need to find the smallest non-negative integer E that gives a valid n (>=3).
    # Let's choose the smallest possible E.
    min_E = 0
    # This gives:
    # n = 4 + 0
    n_for_min_V = 4 + min_E

    print(f"Let's test if a solution exists for V = {min_V}.")
    print(f"To do this, we can choose the minimum possible number of external lines, E = {min_E}.")
    print(f"Plugging V={min_V} and E={min_E} into the equation gives:")
    print(f"n = 4 + E  =>  n = 4 + {min_E} = {n_for_min_V}")
    print(f"This corresponds to a phi^{n_for_min_V} theory, which is a valid interacting scalar field theory.\n")

    print("Therefore, the minimum number of vertices is 1.")
    print("This occurs for a vacuum diagram (E=0) in a phi^4 theory.")
    print("\nThe final equation with the values that yield the minimum is:")
    # Final output showing the numbers in the equation
    print(f"{min_V} = (2 + {min_E}) / ({n_for_min_V} - 2)")

solve_feynman_vertices()
print("\n<<<1>>>")