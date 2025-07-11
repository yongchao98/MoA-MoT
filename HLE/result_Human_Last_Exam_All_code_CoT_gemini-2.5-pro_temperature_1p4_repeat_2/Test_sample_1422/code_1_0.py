def find_minimum_vertices():
    """
    Calculates the minimum number of vertices in a two-loop Feynman diagram
    for a general interacting scalar field theory.
    """
    min_vertices = float('inf')
    best_params = {}

    # L is the number of loops, fixed at 2
    L = 2
    print(f"For a Feynman diagram with L={L} loops, the number of internal lines (I) and vertices (V) are related by:")
    print(f"L = I - V + 1  =>  {L} = I - V + 1  =>  I = V + 1\n")

    print("For a phi^n interaction, the total line ends (n*V) relate to internal (I) and external (E) lines by:")
    print("n * V = 2 * I + E\n")

    print("Combining these gives a formula for V in terms of n and E:")
    print("V = (E + 2) / (n - 2)\n")

    print("Searching for the minimum positive integer V by testing different n and E...")
    print("-" * 30)

    # We search a plausible range for n (interaction type) and E (external lines)
    # n must be >= 3 for an interaction.
    # E is typically even, starting from 0 (vacuum diagrams).
    for n in range(3, 11):  # e.g., phi^3 to phi^10 theory
        for E in range(0, 11, 2):  # e.g., 0 to 10 external lines
            # Numerator and denominator for the formula V = (E + 2) / (n - 2)
            numerator = E + 2
            denominator = n - 2

            if denominator == 0:
                continue

            # Check if V is a positive integer
            if numerator % denominator == 0:
                V = numerator // denominator
                if V > 0 and V < min_vertices:
                    min_vertices = V
                    # Store all sets of parameters that give this minimum
                    best_params = {'V': V, 'n': n, 'E': E, 'num': numerator, 'den': denominator}
                    print(f"Found a potential minimum!")
                    print(f"For a phi^{n} theory (n={n}) with {E} external lines (E={E}):")
                    # Here we output each number in the final equation
                    print(f"V = ({E} + 2) / ({n} - 2) = {best_params['num']} / {best_params['den']} = {best_params['V']}")
                    print("-" * 30)


    print("\n--- Final Result ---")
    if min_vertices == float('inf'):
        print("No solution found in the searched range.")
    else:
        print(f"The minimum number of vertices in a two-loop diagram is {min_vertices}.")
        print("This can be achieved, for example, in a phi^4 theory with 0 external lines (a 'figure-eight' vacuum diagram).")


if __name__ == '__main__':
    find_minimum_vertices()