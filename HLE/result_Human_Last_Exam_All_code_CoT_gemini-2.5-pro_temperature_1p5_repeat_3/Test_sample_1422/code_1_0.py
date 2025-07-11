def find_min_vertices_feynman():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for a standard interacting scalar field theory.
    """

    loops = 2

    print("To find the minimum number of vertices (V) for a two-loop (L=2) Feynman diagram,")
    print("we use two fundamental topological formulas for vacuum diagrams:")
    print("1. Loop Formula: L = I - V + 1, where I is the number of internal lines.")
    print("2. Vertex-Propagator Relation: 2 * I = n * V, where n is the number of lines per vertex.")
    print("\nFirst, we solve for V in terms of L and n:")
    print("From (2), I = (n * V) / 2.")
    print("Substitute I into (1): L = (n * V) / 2 - V + 1")
    print("L - 1 = V * (n/2 - 1)")
    print("V = (L - 1) / (n/2 - 1) = 2 * (L - 1) / (n - 2)")

    print("\nNow, let's plug in L=2 and test for common scalar interactions (n=3, n=4).")

    # Case 1: phi^4 theory (n=4)
    n_phi4 = 4
    # The formula is V = 2 * (L - 1) / (n - 2)
    V_phi4 = 2 * (loops - 1) / (n_phi4 - 2)
    print(f"\nFor a phi^4 theory (n={n_phi4}):")
    print(f"V = 2 * ({loops} - 1) / ({n_phi4} - 2) = 2 * ({loops - 1}) / ({n_phi4 - 2}) = {int(V_phi4)}")
    
    # Case 2: phi^3 theory (n=3)
    n_phi3 = 3
    # The formula is V = 2 * (L - 1) / (n - 2)
    V_phi3 = 2 * (loops - 1) / (n_phi3 - 2)
    print(f"\nFor a phi^3 theory (n={n_phi3}):")
    print(f"V = 2 * ({loops} - 1) / ({n_phi3} - 2) = 2 * ({loops - 1}) / ({n_phi3 - 2}) = {int(V_phi3)}")

    min_vertices = min(V_phi4, V_phi3)
    print(f"\nComparing the possible values ({int(V_phi4)} for phi^4, {int(V_phi3)} for phi^3),")
    print(f"the minimum number of vertices required is {int(min_vertices)}.")
    print("This corresponds to a 'figure-8' tadpole diagram in phi^4 theory.")

find_min_vertices_feynman()