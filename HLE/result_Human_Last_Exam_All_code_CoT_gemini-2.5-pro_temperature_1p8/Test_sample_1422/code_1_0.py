def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for a phi^4 interacting scalar field theory by printing the step-by-step
    derivation.
    """
    
    # Parameters
    L = 2  # Number of loops
    phi_n = 4  # phi^4 theory, so each vertex has 4 lines

    print("Step 1: State the topological formulas for Feynman diagrams.")
    print("----------------------------------------------------------")
    print("The number of loops (L) is given by: L = I - V + 1")
    print(f"The vertex rule for a phi^{phi_n} theory is: {phi_n}*V = 2*I + E")
    print("Where: V = number of vertices, I = number of internal lines, E = number of external lines.\n")

    print("Step 2: Substitute the known value for the number of loops (L).")
    print("-------------------------------------------------------------")
    print(f"We are given L = {L}. Substituting into the loop formula:")
    print(f"{L} = I - V + 1")
    print("Solving for I, we get: I = V + 1\n")

    print("Step 3: Solve the system of equations.")
    print("-------------------------------------")
    print(f"Now, substitute 'I = V + 1' into the vertex rule '{phi_n}*V = 2*I + E':")
    print(f"{phi_n}*V = 2*(V + 1) + E")
    print(f"{phi_n}*V = 2*V + 2 + E")
    print(f"{phi_n-2}*V = 2 + E")
    print(f"V = (2 + E) / {phi_n-2}\n")

    print("Step 4: Find the minimum number of vertices (V).")
    print("------------------------------------------------")
    print("The equation V = (2 + E) / 2 shows that V depends on E.")
    print("To minimize V, we must use the minimum possible value for E.")
    print("The number of external lines (E) must be a non-negative even integer.")
    print("The minimum possible value for E is 0 (a vacuum diagram).\n")
    
    min_E = 0
    # Equation for V: V = 1 + E/2
    min_V = 1 + min_E / 2

    print("Step 5: Calculate the final answer with minimum E.")
    print("-------------------------------------------------")
    print(f"Using E = {min_E}:")
    print(f"V = (2 + {min_E}) / 2")
    # Present the final equation with numbers filled in.
    print(f"Final calculation: {int(min_V)} = (2 + {min_E}) / 2")
    print("\nThis result corresponds to a 'figure-8' vacuum diagram with 1 vertex,")
    print("2 internal lines, and 0 external lines, which forms 2 loops.")
    print(f"\nThus, the minimum number of vertices is {int(min_V)}.")

solve_feynman_vertices()