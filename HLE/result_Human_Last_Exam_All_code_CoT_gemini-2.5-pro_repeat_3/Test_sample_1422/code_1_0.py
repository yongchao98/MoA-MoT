import math

def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in an interacting scalar field theory.
    """
    # --- Step 1: Define the parameters ---
    # Number of loops for a two-loop diagram
    L = 2
    # Number of lines per vertex for the simplest interacting scalar theory (phi^4)
    n = 4
    # Number of external lines for the simplest diagram (a vacuum diagram)
    E = 0

    print("We want to find the minimum number of vertices (V) for a Feynman diagram.")
    print("The topological relations are:")
    print("1) L = I - V + 1         (Loop formula)")
    print("2) n * V = 2 * I + E     (Vertex-line relation)")
    print("\nWe are given the following parameters:")
    print(f" - Number of loops, L = {L}")
    print(f" - Number of lines per vertex, n = {n} (for a phi^4 theory)")
    print(f" - Number of external lines, E = {E} (for a simple vacuum diagram)")

    # --- Step 2: Set up the equation to solve for V ---
    # From equation (1), we can write I as: I = L + V - 1
    # From equation (2), we can write I as: I = (n*V - E) / 2
    # By equating these two expressions for I, we get:
    # L + V - 1 = (n*V - E) / 2

    print("\nBy combining the two relations, we get the equation to solve for V:")
    print("L + V - 1 = (n * V - E) / 2")

    print("\nSubstituting the given values into this equation:")
    print(f"{L} + V - 1 = ({n} * V - {E}) / 2")

    # --- Step 3: Solve the equation ---
    # 2 + V - 1 = (4*V - 0) / 2
    # 1 + V = 2*V
    # 1 = V
    # We can also solve this programmatically from the general form:
    # 2*(L + V - 1) = n*V - E
    # 2L + 2V - 2 = nV - E
    # 2L - 2 + E = nV - 2V
    # 2L - 2 + E = (n - 2)V
    # V = (2L - 2 + E) / (n - 2)
    V = (2 * L - 2 + E) / (n - 2)

    # Check if the result is an integer
    if V.is_integer():
        V = int(V)
        # We can also find the number of internal lines (I) for this diagram
        I = L + V - 1
        print("\nSolving the final equation gives the minimum number of vertices:")
        print(f"V = {V}")
        print(f"\nThis corresponds to a 'figure-eight' diagram with {V} vertex and {I} internal lines.")
    else:
        print(f"\nNo integer solution exists for V with these parameters.")


solve_feynman_vertices()