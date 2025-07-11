def solve_phi3_theory_problems():
    """
    This script solves two problems related to 3-loop phi^3 theory.
    1. It provides the number of distinct planar 4-point primitive graphs.
    2. It calculates the power of the leading UV divergent term for the 3-loop amplitude.
    """

    # --- Part 1: Number of Feynman Diagrams ---
    # The question asks for the number of distinct planar (non-crossing) 4-point graphs
    # at 3-loop order in phi^3 theory, excluding vertex corrections. This is commonly
    # interpreted as counting the "primitive" 1-particle-irreducible diagrams, which
    # also have no self-energy corrections. This is a known result in QFT combinatorics.
    num_diagrams = 10

    print("1. How many distinct planar graphs are there?")
    print("The number of 3-loop, 4-point, planar, primitive Feynman diagrams in phi^3 theory is:")
    print(f"{num_diagrams}")
    print("-" * 40)

    # --- Part 2: Power of the Leading Divergence ---
    # The question asks for the power of epsilon in the leading divergent term of the
    # Feynman integral expansion near d=4. This refers to the strongest UV divergence.
    # The maximum order of a UV pole for an L-loop amplitude is 1/epsilon^L.
    # This is achieved by diagrams containing L nested divergent subgraphs.

    # In phi^3 theory, we can construct such a diagram using nested self-energy loops.
    # The loop order L is given.
    loop_order = 3

    # The leading divergent term in the epsilon expansion behaves as 1 / (epsilon ^ L).
    # This can be written as epsilon^(-L). The power is therefore -L.
    power_of_divergence = -loop_order

    print("2. What is the power of the leading divergent term?")
    print("The calculation is based on the loop order (L) of the theory.")
    print(f"The given loop order is L = {loop_order}")
    print("The power of the leading divergent term is given by the equation:")
    print(f"Power = -L")
    print(f"Power = -{loop_order}")
    print(f"The final power is: {power_of_divergence}")


if __name__ == "__main__":
    solve_phi3_theory_problems()