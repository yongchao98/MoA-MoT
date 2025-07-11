def find_minimum_vertices():
    """
    Calculates and explains the minimum number of vertices in a two-loop
    Feynman diagram for a scalar field theory.
    """
    # 1. Explain the fundamental formulas
    print("To determine the minimum number of vertices in a two-loop Feynman diagram, we use topological relations from graph theory.")
    print("The number of loops (L) in a connected diagram is given by:")
    print("  L = I - V + 1")
    print("where I is the number of internal lines and V is the number of vertices.")

    print("\nFor a theory with an n-point interaction vertex, the lines and vertices are related by:")
    print("  n * V = 2 * I + E")
    print("where E is the number of external lines.")

    # 2. Combine the formulas to solve for V
    print("\nBy combining these formulas, we can solve for V:")
    print("  V = (2*L - 2 + E) / (n - 2)")

    # 3. Apply problem constraints
    L = 2
    # To minimize V, we use the minimum possible number of external lines, E=0 (a vacuum diagram).
    E = 0
    print(f"\nThe problem specifies a two-loop diagram, so L = {L}.")
    print(f"To find the minimum V, we consider diagrams with the minimum number of external lines, so E = {E}.")
    print(f"Plugging these values into the formula for V:")
    
    numerator = 2 * L - 2 + E
    print(f"  V = (2*{L} - 2 + {E}) / (n - 2) = {numerator} / (n - 2)")

    # 4. Find the minimum integer solution for V by checking different interaction types (n)
    print("\nNow we check different interaction types (values of n) to find the minimum positive integer for V.")

    # Case n=3 (e.g., phi^3 theory)
    n3 = 3
    denominator3 = n3 - 2
    V3 = numerator / denominator3
    print(f"\nFor a cubic interaction (n = {n3}):")
    print(f"  V = {numerator} / ({n3} - 2) = {numerator} / {denominator3} = {int(V3)}")
    
    # Case n=4 (e.g., phi^4 theory)
    n4 = 4
    denominator4 = n4 - 2
    V4 = numerator / denominator4
    print(f"\nFor a quartic interaction (n = {n4}):")
    print(f"  V = {numerator} / ({n4} - 2) = {numerator} / {denominator4} = {int(V4)}")

    print("\nFor any interaction with n > 4, the denominator (n-2) would be greater than 2, making V a fraction less than 1. Since V must be a positive integer, these cases are not possible for a vacuum diagram.")

    final_answer = int(V4)
    print(f"\nTherefore, the minimum number of vertices in a two-loop diagram for an interacting scalar field theory is {final_answer}.")

# Call the function to execute the logic
find_minimum_vertices()