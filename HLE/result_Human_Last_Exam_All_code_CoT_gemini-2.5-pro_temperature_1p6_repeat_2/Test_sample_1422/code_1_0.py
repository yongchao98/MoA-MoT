import math

def calculate_min_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in an interacting scalar field theory.
    """
    #
    # We start with two fundamental topological relations for a connected Feynman diagram:
    # 1) L = I - V + 1
    #    L = number of loops
    #    I = number of internal lines (propagators)
    #    V = number of vertices (interactions)
    #
    # 2) n * V = 2 * I + E
    #    n = number of fields meeting at each vertex (e.g., n=4 for a phi^4 theory)
    #    E = number of external lines (incoming/outgoing particles)
    #
    # Our goal is to solve for V.

    # Step 1: Define the knowns for our specific problem.
    # We want a two-loop diagram.
    L = 2
    # To find the absolute minimum number of vertices, we consider the simplest
    # possible process, which is a vacuum diagram with no external particles.
    E = 0
    # Let's consider the common phi^4 interaction, which often allows for the
    # most compact diagrams. For a phi^4 interaction, n=4.
    n = 4

    # Step 2: Combine the equations to solve for V.
    # From equation (1), we can express I in terms of L and V:
    # I = L + V - 1
    #
    # Now, substitute this expression for I into equation (2):
    # n * V = 2 * (L + V - 1) + E
    # n * V = 2*L + 2*V - 2 + E
    #
    # Now, group the V terms:
    # n*V - 2*V = 2*L - 2 + E
    # V * (n - 2) = 2*L - 2 + E
    #
    # Finally, isolate V:
    # V = (2*L - 2 + E) / (n - 2)
    # This is the general formula we will use.

    # Step 3: Plug in the values and calculate the result.
    numerator = 2 * L - 2 + E
    denominator = n - 2
    
    # Check for division by zero, which happens for n=2 theories.
    if denominator == 0:
        print(f"The calculation is not possible for n={n}.")
        return

    V = numerator / denominator

    # Let's print the full calculation step-by-step.
    print("The formula for the number of vertices (V) is derived as: V = (2*L - 2 + E) / (n - 2)")
    print("-" * 30)
    print(f"Given values:")
    print(f"Number of loops (L) = {L}")
    print(f"Number of external lines (E) = {E} (minimized)")
    print(f"Number of fields per vertex (n) = {n} (for phi^4 theory)")
    print("-" * 30)
    print("Calculation:")
    print(f"V = (2 * {L} - 2 + {E}) / ({n} - 2)")
    print(f"V = ({2*L} - 2 + {E}) / ({n-2})")
    print(f"V = {numerator} / {denominator}")
    print(f"V = {math.ceil(V)}") # Vertices must be integers
    print("-" * 30)
    print("This result means a two-loop diagram is possible with just one vertex in a phi^4 theory.")
    print("This diagram is often called a 'figure-eight' or 'double-scoop' vacuum diagram.")

calculate_min_vertices()