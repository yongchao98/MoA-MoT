import math

def calculate_min_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in an interacting scalar field theory.
    """
    # We are given the number of loops.
    L = 2

    # To find the minimum number of vertices, we consider the simplest possible diagram.
    # 1. Simplest stable interaction: phi^4 theory, where 4 lines meet at each vertex.
    k = 4
    # 2. Simplest process: a vacuum diagram with no external particles.
    E = 0

    # The two key topological formulas for a Feynman diagram are:
    # 1) L = I - V + 1  (where I is internal lines, V is vertices)
    # 2) k * V = 2 * I + E

    # We can solve these for V. From (1), we get I = L + V - 1.
    # Substitute this into (2):
    # k * V = 2 * (L + V - 1) + E
    # k * V = 2*L + 2*V - 2 + E
    # (k - 2) * V = 2*L - 2 + E
    # V = (2*L - 2 + E) / (k - 2)

    # Now, we plug in our values to find V.
    numerator = 2 * L - 2 + E
    denominator = k - 2
    V = numerator / denominator

    # The number of vertices must be an integer.
    V = int(V)

    # We can also find the number of internal lines for this diagram.
    I = L + V - 1

    # --- Output the explanation and result ---
    print("To find the minimum number of vertices in a two-loop Feynman diagram, we use topological formulas.")
    print("\nLet L=loops, V=vertices, I=internal lines, E=external lines, k=lines per vertex.")
    print("The governing equations are:")
    print("  1) L = I - V + 1")
    print("  2) k * V = 2 * I + E")

    print("\nWe are given L=2. For the *minimum* number of vertices, we assume the simplest case:")
    print(f"  - A phi^4 interaction, so k = {k}")
    print(f"  - A vacuum diagram, so E = {E}")

    print("\nSolving for V, we get the formula: V = (2*L - 2 + E) / (k - 2)")
    print("\nPlugging in the numbers:")
    print(f"  V = (2 * {L} - 2 + {E}) / ({k} - 2)")
    print(f"  V = ({2*L} - 2 + {E}) / ({k-2})")
    print(f"  V = {numerator} / {denominator}")
    print(f"  V = {V}")

    print(f"\nThis minimal diagram has V={V} vertex and I={I} internal lines.")
    print("This corresponds to a 'figure-eight' vacuum diagram, where a single 4-point vertex has its lines connected in pairs to form two loops.")
    print(f"Checking the loop formula: L = I - V + 1 = {I} - {V} + 1 = {I-V+1}. This is correct.")
    print("\nTherefore, the minimum number of vertices is 1.")

if __name__ == '__main__':
    calculate_min_vertices()