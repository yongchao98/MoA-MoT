def solve_feynman_vertices():
    """
    Calculates and explains the minimum number of vertices in a two-loop
    Feynman diagram for an interacting scalar field theory.
    """
    # The number of loops is given.
    L = 2

    print("To determine the minimum number of vertices, we use topological relations from graph theory applied to Feynman diagrams.")
    
    print("\nStep 1: Relate loops (L), internal lines (I), and vertices (V).")
    print("The formula is: L = I - V + 1")
    print(f"For a two-loop diagram, we have L = {L}. Substituting this in:")
    print(f"  {L} = I - V + 1")
    print("  This gives us a relation between I and V: I = V + 1\n")

    print("Step 2: Relate vertices and lines for a specific theory.")
    print("For an interacting scalar field theory where the interaction vertex involves 'n' fields (a 'phi-n' theory), the total number of line-ends (n*V) must equal twice the number of internal lines plus the number of external lines (E).")
    print("The formula is: n * V = 2 * I + E\n")
    
    print("Step 3: Combine the formulas and solve for V.")
    print("We substitute 'I = V + 1' into the second equation:")
    print("  n * V = 2 * (V + 1) + E")
    print("  n * V = 2*V + 2 + E")
    print("  V * (n - 2) = 2 + E")
    print("  V = (2 + E) / (n - 2)\n")

    print("Step 4: Find the minimum integer V by checking simple theories.")
    print("To find the minimum V, we should use the simplest diagrams (minimum external lines, E) and simplest interactions (minimum n).")
    print("An interaction requires n >= 3. The minimum number of external lines is E = 0 (a vacuum diagram).\n")
    
    # --- Case 1: phi^3 Theory ---
    n1 = 3
    e1 = 0
    print(f"First, consider the phi^{n1} theory (n={n1}) with E={e1} external lines:")
    v1_num = 2 + e1
    v1_den = n1 - 2
    v1 = v1_num / v1_den
    print(f"  The equation for V is: V = ({2} + {e1}) / ({n1} - {2})")
    print(f"  Which evaluates to: V = {v1_num} / {v1_den} = {int(v1)}")
    print(f"  So, for a phi^3 theory, the minimum number of vertices is 2.\n")

    # --- Case 2: phi^4 Theory ---
    n2 = 4
    e2 = 0
    print(f"Next, consider the phi^{n2} theory (n={n2}) with E={e2} external lines:")
    v2_num = 2 + e2
    v2_den = n2 - 2
    v2 = v2_num / v2_den
    print(f"  The equation for V is: V = ({2} + {e2}) / ({n2} - {2})")
    print(f"  Which evaluates to: V = {v2_num} / {v2_den} = {int(v2)}")
    print(f"  For a phi^4 theory, a diagram with just 1 vertex is possible.\n")

    # --- Conclusion ---
    min_v = min(v1, v2)
    print("Comparing the simplest interaction theories (phi^3 and phi^4), the absolute minimum possible number of vertices is the smaller result.")
    print("This occurs for a phi^4 vacuum diagram, which looks like a figure-eight.")
    print("\nFinal Answer:")
    print(f"The minimum number of vertices is {int(min_v)}.")


solve_feynman_vertices()