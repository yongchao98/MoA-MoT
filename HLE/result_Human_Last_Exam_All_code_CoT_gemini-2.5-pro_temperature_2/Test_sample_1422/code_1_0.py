import math

def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram.
    """
    print("### Finding the Minimum Number of Vertices for a Two-Loop Diagram ###")
    print("-" * 60)

    # Step 1 & 2: Use the loop formula
    print("Step 1: Use the topological formula for Feynman diagrams.")
    print("The number of loops (L) is given by L = I - V + 1, where:")
    print("  I = number of internal lines")
    print("  V = number of vertices")
    
    L = 2
    print(f"\nFor a two-loop diagram, L = {L}. So, the equation becomes:")
    print(f"{L} = I - V + 1")
    print("Rearranging this, we find the relationship: I = V + 1\n")
    print("-" * 60)
    
    # Step 3: Use the vertex formula
    print("Step 2: Relate vertices to lines.")
    print("For a theory with only one type of n-point vertex, the total number of line 'ends' is n * V.")
    print("These ends connect to internal lines (twice each) or external lines (once each).")
    print("This gives the formula: n * V = 2 * I + E, where:")
    print("  n = number of lines at each vertex")
    print("  E = number of external lines (incoming/outgoing particles)\n")
    print("-" * 60)
    
    # Step 4: Combine and analyze
    print("Step 3: Combine the formulas to solve for V.")
    print("Substitute I = V + 1 into the second formula:")
    print("n * V = 2 * (V + 1) + E")
    print("n * V = 2*V + 2 + E")
    print("(n - 2) * V = 2 + E")
    print("So, V = (2 + E) / (n - 2)\n")
    print("To minimize V, we must minimize E (number of external lines) and use an interaction n > 2.")
    print("The minimal number of external lines is E = 0 (a 'vacuum bubble' diagram).")
    print("-" * 60)

    # Step 5 & 6: Test specific theories
    print("Step 4: Calculate V for common interaction theories with E = 0.\n")
    
    # Case 1: phi^3 theory
    n_3 = 3
    E_min = 0
    print(f"--- Case: phi^3 Theory (n = {n_3}) ---")
    numerator_3 = 2 + E_min
    denominator_3 = n_3 - 2
    v_3 = numerator_3 / denominator_3
    print(f"V = (2 + E) / (n - 2) = ({2} + {E_min}) / ({n_3} - {denominator_3}) = {int(v_3)}")
    print(f"The minimum vertices for a phi^3 theory is {int(v_3)}.\n")

    # Case 2: phi^4 theory
    n_4 = 4
    print(f"--- Case: phi^4 Theory (n = {n_4}) ---")
    print(f"(Note: For a phi^4 vertex, E must be even. E={E_min} is valid.)")
    numerator_4 = 2 + E_min
    denominator_4 = n_4 - 2
    v_4 = numerator_4 / denominator_4
    print(f"V = (2 + E) / (n - 2) = ({2} + {E_min}) / ({n_4} - {denominator_4}) = {int(v_4)}")
    print(f"The minimum vertices for a phi^4 theory is {int(v_4)}.\n")

    # Final conclusion
    print("-" * 60)
    print("CONCLUSION:")
    print("Comparing the cases, the minimum number of vertices is achieved with the phi^4 theory.")
    final_answer = int(min(v_3, v_4))
    print(f"The absolute minimum number of vertices required is {final_answer}.")
    print("This corresponds to a phi^4 vacuum diagram with 1 vertex and 2 internal lines in a 'figure-8' configuration.")

solve_feynman_vertices()