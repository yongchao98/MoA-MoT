def calculate_min_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in a scalar phi-4 theory.
    """
    # Number of loops specified in the problem
    num_loops = 2

    # Explain the context and formulas
    print("To find the minimum number of vertices, we analyze an interacting scalar field theory.")
    print("The phi-4 interaction (4 lines per vertex) typically provides the minimum vertex count.")
    print("We use two key formulas:")
    print("1. Topological formula: L = I - V + 1")
    print("2. Phi-4 vertex rule: 4V = 2I + E")
    print("Where L=loops, V=vertices, I=internal lines, E=external lines.")
    print("\nBy combining these formulas, we can derive a direct relationship between V, L, and E:")
    print("V = L - 1 + E / 2")

    # To minimize vertices (V) for a fixed number of loops (L), we must use the
    # minimum possible number of external lines (E).
    # A diagram for a physical process (like self-energy) needs at least 2 external lines.
    min_external_lines = 2

    print(f"\nThe problem specifies a two-loop diagram, so L = {num_loops}.")
    print(f"To minimize V, we use the minimum number of external lines, E = {min_external_lines}.")
    
    # Perform the calculation
    min_vertices = num_loops - 1 + min_external_lines / 2

    # Output the final calculation and result
    print("\n--- Calculation ---")
    print(f"Minimum Vertices (V) = {num_loops} - 1 + {min_external_lines} / 2")
    # To show the steps clearly
    step_1 = num_loops - 1
    step_2 = min_external_lines / 2
    print(f"V = {step_1} + {step_2}")
    print(f"V = {int(min_vertices)}")
    print("-------------------")
    print(f"\nThe minimum number of vertices is {int(min_vertices)}.")

calculate_min_vertices()