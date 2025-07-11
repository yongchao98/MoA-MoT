def solve_feynman_diagram_problems():
    """
    This function solves the two specified problems regarding 3-loop Feynman diagrams
    in phi^3 theory.
    """

    # --- Problem 1: Number of distinct planar graphs ---
    # As explained in the reasoning, the number of distinct graphs can be determined
    # by enumerating the ways to insert two self-energy bubbles into a 1-loop box graph,
    # which yields 3 non-isomorphic graphs.
    num_graphs = 3
    print(f"1. The number of distinct planar graphs is: {num_graphs}")
    print("-" * 20)

    # --- Problem 2: Power of the leading divergent term ---
    # For a massless on-shell scattering amplitude at L-loop order, the leading
    # infrared divergence in d=4-2*epsilon dimensions has a pole of order 2*L.
    num_loops = 3
    leading_divergence_power = 2 * num_loops

    print("2. The power of the leading divergent term of the epsilon expansion is calculated as follows:")
    # The final print statement is formatted to show the equation as requested.
    print(f"Power = 2 * L")
    print(f"Power = 2 * {num_loops} = {leading_divergence_power}")


if __name__ == "__main__":
    solve_feynman_diagram_problems()
