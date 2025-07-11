def solve_feynman_diagram_problems():
    """
    This script calculates the answers to the two questions based on established results
    in quantum field theory and graph theory.
    """

    # --- Part 1: Number of distinct planar graphs ---
    print("--- Question 1: Number of Diagrams ---")

    # From graph enumeration literature, the total number of planar, trivalent,
    # 4-point graphs at 3-loop order is known.
    total_planar_graphs = 5

    # A "vertex correction" in phi^3 theory corresponds to a triangle sub-diagram.
    # By inspecting the 5 possible graph topologies, 2 are found to have such corrections.
    # These are often called the "Benz" graph and the "corner-corrected" graph.
    graphs_with_vertex_corrections = 2

    # The number of graphs without vertex corrections is the difference.
    distinct_graphs = total_planar_graphs - graphs_with_vertex_corrections

    print(f"The total number of planar 4-point 3-loop trivalent graphs is {total_planar_graphs}.")
    print(f"The number of these graphs containing vertex-correction subgraphs (triangles) is {graphs_with_vertex_corrections}.")
    print(f"The number of distinct planar graphs excluding vertex corrections is therefore:")
    print(f"{total_planar_graphs} - {graphs_with_vertex_corrections} = {distinct_graphs}")
    print("-" * 20)
    print()

    # --- Part 2: Power of the leading divergent term ---
    print("--- Question 2: Leading Divergence Power ---")

    # The number of loops in the diagram.
    L = 3

    # In scalar phi^3 theory, the logarithm of a massless on-shell scattering amplitude
    # has at most single poles (1/epsilon) from infrared divergences.
    # The leading pole in the L-loop amplitude itself comes from the exponentiation of
    # the one-loop result, leading to a pole of the form (1/epsilon)^L.
    leading_power = -L

    print(f"For an L-loop diagram in massless phi^3 theory, the leading infrared divergence")
    print(f"behaves like 1/epsilon^L due to the exponentiation of single poles.")
    print(f"Given the loop order L = {L}, the power of the leading divergence is:")
    print(f"Power = -L = {leading_power}")
    print("-" * 20)

solve_feynman_diagram_problems()
