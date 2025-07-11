def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in an interacting scalar field theory.
    """
    # 1. State the problem and the assumptions for finding the minimum.
    print("Problem: Find the minimum number of vertices (V) in a two-loop Feynman diagram.")
    print("-" * 75)
    print("To find the minimum, we will analyze the simplest possible configuration:")
    print(" - Theory: An interacting scalar phi^4 theory, where each vertex has n=4 lines.")
    print(" - Process: A vacuum diagram, which has E=0 external lines.")
    print(" - Loops: We are looking for a two-loop diagram, so L=2.")
    print("\n")

    # 2. State the general formulas.
    print("The key topological relations for a connected Feynman diagram are:")
    print(" (1)  L = I - V + 1")
    print(" (2)  n * V = 2 * I + E")
    print("      where L=loops, I=internal lines, V=vertices, n=lines/vertex, E=external lines.")
    print("\n")

    # 3. Substitute values and derive two expressions for I.
    L = 2
    n = 4
    E = 0

    print("Substituting our values (L=2, n=4, E=0) into the formulas:")

    print("\nFrom formula (1):")
    print(f"  {L} = I - V + 1")
    print(f"  I = V + {L} - 1")
    print(f"  I = V + {L - 1}")

    print("\nFrom formula (2):")
    print(f"  {n} * V = 2 * I + {E}")
    print(f"  {n} * V = 2 * I")
    print(f"  I = ({n}/2) * V")
    print(f"  I = {int(n/2)} * V")
    print("\n")

    # 4. Equate the expressions and solve for V.
    print("Now we set the two expressions for I equal to each other to solve for V:")
    # The final equation is V + 1 = 2 * V
    print("Equation: V + 1 = 2 * V")
    print("\nSolving for V:")
    print("  1 = 2 * V - V")
    final_V = 1
    print(f"  1 = {final_V} * V")
    print(f"  V = {final_V}")
    print("-" * 75)

    print(f"\nThe minimum number of vertices required is {final_V}.")
    print("This corresponds to a 'figure-eight' vacuum diagram, where a single 4-point vertex")
    print("has its lines connected to form two loops.")

# Run the calculation and print the result.
solve_feynman_vertices()