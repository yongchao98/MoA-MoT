def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for an interacting scalar field theory (phi^4).
    """

    print("To determine the minimum number of vertices, we will use the topological rules for Feynman diagrams, assuming a phi^4 interaction theory.")
    print("The variables are: V = number of vertices, I = number of internal lines, E = number of external lines, and L = number of loops.")
    print("-" * 60)

    # Step 1: State the general formula for the number of loops.
    print("1. The number of loops (L) in any connected diagram is given by:")
    print("   L = I - V + 1")
    print("\nFor a two-loop diagram, we set L = 2:")
    print("   2 = I - V + 1")
    print("   Rearranging this gives us a relation between I and V:")
    print("   I = V + 1")
    print("-" * 60)

    # Step 2: State the formula for vertices in phi^4 theory.
    print("2. For a phi^4 interaction, each vertex has four lines attached. The total number")
    print("   of line 'ends' (4*V) must equal the sum of ends forming internal lines (2*I)")
    print("   and external lines (E). This gives the relation:")
    print("   4 * V = 2 * I + E")
    print("-" * 60)

    # Step 3: Solve the system of equations.
    print("3. Now, we substitute the expression for I from step 1 into the equation from step 2:")
    print("   4 * V = 2 * (V + 1) + E")
    print("   4 * V = 2 * V + 2 + E")
    print("\n   Subtracting 2*V from both sides gives:")
    print("   2 * V = 2 + E")
    print("\n   Finally, solving for V:")
    print("   V = 1 + E / 2")
    print("-" * 60)

    # Step 4: Find the minimum value.
    print("4. The equation V = 1 + E/2 shows that to minimize the number of vertices (V),")
    print("   we must minimize the number of external lines (E). The smallest possible")
    print("   number of external lines for a diagram is E = 0 (a vacuum diagram).")
    print("-" * 60)

    # Step 5: Calculate the final answer.
    E_min = 0
    V_min_val_part1 = 1
    V_min_val_part2_num = E_min
    V_min_val_part2_den = 2
    V_min = V_min_val_part1 + V_min_val_part2_num / V_min_val_part2_den

    print("5. Plugging E = 0 into our final equation gives the minimum number of vertices:")
    print(f"   V_min = {V_min_val_part1} + {V_min_val_part2_num} / {V_min_val_part2_den}")
    print(f"   V_min = {int(V_min)}")
    print("-" * 60)

    print(f"Therefore, the minimum number of vertices in a two-loop Feynman diagram for this theory is {int(V_min)}.")

solve_feynman_vertices()