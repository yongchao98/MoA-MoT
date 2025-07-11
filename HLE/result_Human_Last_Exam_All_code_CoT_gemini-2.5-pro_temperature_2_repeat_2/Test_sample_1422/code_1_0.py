def solve_feynman_diagram_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for a scalar field theory.
    """

    # We are given the number of loops in the diagram.
    L = 2
    print(f"Goal: Find the minimum vertices (V) for a diagram with {L} loops.\n")

    # The standard model for a simple interacting scalar field is the phi-4 theory,
    # where each vertex has 4 lines connected to it.
    n = 4
    print(f"We will analyze the most common case: φ^{n} theory, where n = {n}.")
    print("-" * 50)

    # We use two fundamental relationships for Feynman diagrams.
    print("Step 1: Use the topological formula relating Loops (L), Internal Lines (I), and Vertices (V).")
    print("Formula: L = I - V + 1")
    print(f"For L = {L}, we have: {L} = I - V + 1")
    print("Rearranging for I, we get: I = V + 1\n")

    print("Step 2: Use the vertex rule relating n, V, I, and External Lines (E).")
    print("Formula: n * V = 2 * I + E")
    print(f"Substituting n = {n}, this becomes: {n} * V = 2 * I + E\n")

    print("Step 3: Combine the two formulas.")
    print("Substitute 'I = V + 1' into the vertex rule:")
    print(f"{n} * V = 2 * (V + 1) + E")
    print(f"{n} * V = 2 * V + 2 + E")
    print(f"({n} - 2) * V = 2 + E")
    
    # Substituting the value of n
    print(f"({n-2}) * V = 2 + E\n")

    print("Step 4: Solve for the minimum V.")
    print("To find the minimum number of vertices (V), we must choose the number")
    print("of external lines (E) that makes V as small as possible.")
    print("The minimum possible number of external lines is E = 0 (a vacuum diagram).\n")
    
    E = 0
    print(f"Plugging in E = {E}:")
    # Equation is (n-2)*V = 2 + E
    lhs_coeff = n - 2
    rhs = 2 + E
    print(f"{lhs_coeff} * V = {rhs}")
    
    min_V = rhs // lhs_coeff
    print(f"V = {rhs} / {lhs_coeff}")
    print(f"V = {min_V}\n")

    print("-" * 50)
    print(f"Conclusion: The minimum number of vertices required for a two-loop diagram is {min_V}.")
    print("This corresponds to the 'figure-eight' vacuum diagram in φ^4 theory.")

solve_feynman_diagram_vertices()