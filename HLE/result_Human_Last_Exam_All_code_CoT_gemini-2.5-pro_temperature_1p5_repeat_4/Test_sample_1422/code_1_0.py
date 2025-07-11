def calculate_minimum_vertices():
    """
    Calculates the minimum number of vertices for a two-loop diagram
    in an interacting scalar field theory.
    """
    # Number of loops as given by the problem.
    L = 2

    # For the interaction type, we choose k=4, which corresponds to the common
    # and simple phi^4 theory. At each vertex, k lines meet.
    k = 4

    # To minimize the number of vertices V, we must choose the minimum
    # possible number of external lines E. The absolute minimum is E=0,
    # which represents a vacuum diagram.
    E = 0

    # The number of vertices V can be derived from the topological relations of
    # Feynman diagrams. The general formula is V = (2*L - 2 + E) / (k - 2).
    # For L=2, this simplifies to V = (2 + E) / (k - 2).
    # We use integer division // as the result must be an integer.
    V = (2 + E) // (k - 2)

    # Output the explanation and the final calculation.
    print("To find the minimum number of vertices (V) in a two-loop (L=2) diagram:")
    print("We use the general formula: V = (2 * L - 2 + E) / (k - 2)")
    print(f"For L={L}, the formula becomes: V = (2 + E) / (k - 2)\n")

    print("We can achieve the minimum V by choosing:")
    print(f"- A simple interaction type, like phi^k with k = {k}")
    print(f"- The minimum number of external lines, E = {E}\n")

    print("Now, we substitute these values into the equation:")
    print(f"V = (2 + {E}) / ({k} - 2)")
    print(f"V = {2 + E} / {k - 2}")
    print(f"V = {V}")

    print(f"\nThe minimum number of vertices in a two-loop Feynman diagram for an interacting scalar field theory is {V}.")

calculate_minimum_vertices()