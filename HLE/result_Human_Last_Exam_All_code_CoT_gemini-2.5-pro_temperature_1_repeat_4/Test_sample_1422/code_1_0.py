def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for a phi^4 interacting scalar field theory.
    """
    # Number of loops required
    L = 2

    print("We want to find the minimum number of vertices (V) for a two-loop Feynman diagram.")
    print("For a phi^4 theory, we use two fundamental equations:")
    print("1. Loop formula: L = I - V + 1")
    print("2. Vertex rule: 4V = 2I + E")
    print("Where L is loops, I is internal lines, V is vertices, and E is external lines.\n")

    print("First, we combine the equations to solve for V:")
    print("From (1), we express I: I = L + V - 1")
    print("Substitute this into (2): 4V = 2 * (L + V - 1) + E")
    print("Solving for V gives: V = L - 1 + E / 2\n")

    print(f"We are given L = {L}.")
    print("To find the minimum V, we must use the minimum possible number of external lines (E).")
    print("The minimum is E = 0 (a vacuum diagram).\n")

    # Minimum number of external lines
    E = 0

    # Calculate the minimum number of vertices
    min_V = L - 1 + E / 2

    # The result must be an integer
    min_V = int(min_V)

    print("--- Final Calculation ---")
    print(f"V = L - 1 + E / 2")
    print(f"V = {L} - 1 + {E} / 2")
    print(f"V = {min_V}")
    print("-------------------------\n")
    print(f"The minimum number of vertices is {min_V}.")

solve_feynman_vertices()
<<<1>>>