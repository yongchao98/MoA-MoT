def find_minimum_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for a phi^4 interacting scalar field theory.
    """
    print("This program finds the minimum number of vertices for a two-loop Feynman diagram.")
    print("Assumption: We are using an interacting scalar field theory with a 4-point vertex (phi^4 theory).\n")

    print("Step 1: Define the topological relations.")
    print("Let V = number of vertices, I = number of internal lines, E = number of external lines, L = number of loops.")
    print("Formula 1 (Loops): L = I - V + 1")
    print("Formula 2 (Vertices for phi^4 theory): 4 * V = 2 * I + E\n")

    print("Step 2: Solve for V given L=2.")
    print("We are given L = 2. From Formula 1, we get:")
    print("2 = I - V + 1  =>  I = V + 1\n")

    print("Now, substitute I = V + 1 into Formula 2:")
    print("4*V = 2*(V + 1) + E")
    print("4*V = 2*V + 2 + E")
    print("2*V = 2 + E")
    print("V = 1 + E / 2\n")

    print("Step 3: Minimize V by minimizing E.")
    print("The number of external lines, E, must be a non-negative even integer.")
    print("To find the minimum number of vertices V, we must use the minimum possible value for E.")
    min_E = 0
    print(f"The minimum possible value for E is {min_E} (corresponding to a vacuum diagram).\n")
    
    print("Step 4: Calculate the final minimum number of vertices.")
    min_V = 1 + min_E / 2

    # The user requested to output each number in the final equation.
    print("The final calculation is:")
    print(f"V_min = 1 + {min_E} / 2 = {int(min_V)}")

    # Store the final result
    result = int(min_V)
    print(f"\nThus, the minimum number of vertices in a two-loop Feynman diagram for this theory is {result}.")

    # Output the final answer in the requested format
    print(f"\n<<<{result}>>>")

if __name__ == '__main__':
    find_minimum_vertices()
