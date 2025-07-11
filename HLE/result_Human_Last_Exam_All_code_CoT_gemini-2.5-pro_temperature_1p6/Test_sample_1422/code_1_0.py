def calculate_min_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in an interacting scalar field theory by considering common interaction types.
    """
    print("The number of vertices (V) in a Feynman diagram is related to the number of loops (L),")
    print("external lines (E), and the type of interaction vertex (n-point) by the formula:")
    print("V * (n - 2) = 2*L - 2 + E\n")

    # Given value for loops
    L = 2
    print(f"We are given L = {L}, so the formula becomes: V * (n - 2) = 2 + E.\n")

    # --- Case 1: phi^4 theory ---
    n_phi4 = 4
    # To minimize V, we take the minimum possible number of external lines, E=0 (a vacuum diagram).
    E_min = 0
    
    print(f"First, let's consider a phi^4 theory, where each vertex connects 4 lines (n = {n_phi4}).")
    print(f"To minimize V, we should use the minimum number of external lines, so we choose E = {E_min}.")
    print("Plugging the values into the formula:")
    
    # Calculate each part of the equation to show the steps
    V_mult_phi4 = n_phi4 - 2
    rhs_phi4 = 2 + E_min
    V_phi4 = rhs_phi4 / V_mult_phi4

    print(f"V * ({n_phi4} - 2) = 2 + {E_min}")
    print(f"V * {V_mult_phi4} = {rhs_phi4}")
    print(f"V = {rhs_phi4} / {V_mult_phi4}")
    print(f"V = {int(V_phi4)}\n")


    # --- Case 2: phi^3 theory ---
    n_phi3 = 3
    # Again, we take E=0 to minimize V.
    
    print(f"Next, let's consider a phi^3 theory, where each vertex connects 3 lines (n = {n_phi3}).")
    print(f"Again, we choose E = {E_min} to find the minimum V.")
    print("Plugging the values into the formula:")

    # Calculate each part of the equation
    V_mult_phi3 = n_phi3 - 2
    rhs_phi3 = 2 + E_min
    V_phi3 = rhs_phi3 / V_mult_phi3

    print(f"V * ({n_phi3} - 2) = 2 + {E_min}")
    print(f"V * {V_mult_phi3} = {rhs_phi3}")
    print(f"V = {rhs_phi3} / {V_mult_phi3}")
    print(f"V = {int(V_phi3)}\n")

    # --- Conclusion ---
    final_answer = min(int(V_phi4), int(V_phi3))
    print(f"Comparing the results ({int(V_phi4)} for phi^4 vs. {int(V_phi3)} for phi^3), the minimum possible number of vertices is {final_answer}.")
    print("This corresponds to a 'figure-eight' vacuum diagram in a phi^4 theory.")

calculate_min_vertices()