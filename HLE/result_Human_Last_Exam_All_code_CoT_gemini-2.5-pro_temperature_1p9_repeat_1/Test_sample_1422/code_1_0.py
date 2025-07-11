def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for a phi-4 interacting scalar field theory.
    """

    # The number of loops (L) is given.
    L = 2
    
    print("This program calculates the minimum number of vertices (V) for a two-loop Feynman diagram.")
    print("We assume a phi-4 scalar field theory, where each vertex has 4 lines.")
    print("-" * 50)
    
    # Explain the two fundamental equations
    print("We start with two key equations:")
    print("1. Topological formula for loops: L = I - V + 1")
    print("2. Vertex connection rule for phi-4 theory: 4*V = 2*I + E")
    print("   (I = internal lines, E = external lines)")
    print("-" * 50)
    
    # Step 1: Solve the first equation for I
    print(f"Step 1: Using L = {L}, we rearrange the first equation to solve for I:")
    print(f"   {L} = I - V + 1")
    print("   => I = V + 1")
    print("-" * 50)
    
    # Step 2: Substitute I into the second equation
    print("Step 2: Substitute I = V + 1 into the second equation:")
    print("   4*V = 2*(V + 1) + E")
    print("   4*V = 2*V + 2 + E")
    print("   2*V = 2 + E")
    print("   V = 1 + E / 2")
    print("-" * 50)

    # Step 3: Find the minimum V by minimizing E
    print("Step 3: To find the minimum V, we must use the minimum possible value for E.")
    print("The number of external lines E can be 0 (for a vacuum diagram).")
    
    min_E = 0
    print(f"We set E = {min_E} to find the minimum number of vertices.")
    print("-" * 50)
    
    # Final Calculation
    print("Final Calculation:")
    # The final equation requires printing each number
    final_v = 1 + min_E / 2
    print(f"V = 1 + {min_E} / 2")
    print(f"V = {int(final_v)}")
    
    print("\nConclusion:")
    print(f"The minimum number of vertices is {int(final_v)}. This corresponds to a 'figure-eight' vacuum diagram.")

solve_feynman_vertices()