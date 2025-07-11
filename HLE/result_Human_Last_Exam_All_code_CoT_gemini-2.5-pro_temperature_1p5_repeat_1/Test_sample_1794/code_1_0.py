import math

def solve_phi3_theory_problems():
    """
    This script provides solutions to two questions about 3-loop, 4-point
    scalar phi^3 theory. It explains the reasoning and prints the numerical answers.
    """

    # --- Part 1: Counting Planar Graphs ---
    print("--- Question 1: Number of distinct planar graphs ---")
    
    # In phi^3 theory, the number of planar, 4-point Feynman diagrams for a fixed
    # cyclic ordering of external legs is given by a known combinatorial sequence.
    # For L=1, 2, 3, and 4 loops, the number of diagrams are 1, 2, 5, and 14, respectively.
    # For this problem at L=3 loops, the total number of planar topologies is 5.
    
    # The question asks to exclude diagrams with "vertex corrections". A vertex
    # correction diagram is formed by replacing a vertex in a lower-loop diagram
    # with a 1-loop triangle. The 5 standard planar topologies are either primitive
    # (the "triple-ladder" graph) or contain self-energy insertions (a loop on a
    # propagator). These structures are distinct from vertex correction structures.
    # Therefore, we count all 5 of these diagrams.

    num_graphs = 5
    
    print("The number of 3-loop, 4-point planar graphs, excluding vertex corrections, is:")
    print(num_graphs)
    print("-" * 20)


    # --- Part 2: Power of the Leading Divergent Term ---
    print("--- Question 2: Power of the leading divergent term ---")

    # The diagram is for massless on-shell particles in d = 4 - 2*epsilon dimensions.
    # The leading divergence (highest pole in epsilon) can be from Ultraviolet (UV)
    # or Infrared (IR) regions of the loop momentum integration.

    # Maximum UV pole order:
    # In 4D phi^3 theory, 4-point functions are primitively UV finite. Divergences
    # arise only from self-energy subgraphs. A 3-loop 1-Particle-Irreducible (1PI)
    # diagram can have a maximum UV pole of 1/epsilon^2. This comes from diagrams
    # containing two non-overlapping 1-loop self-energy subgraphs, or one nested
    # 2-loop self-energy subgraph.
    max_uv_pole_order = 2

    # Maximum IR pole order:
    # For massless on-shell scattering, IR divergences arise. A general result for
    # scalar field theories is that an L-loop amplitude can have IR poles of at most
    # order L. For L=3 loops, the maximum IR pole is 1/epsilon^3.
    max_ir_pole_order = 3
    
    # The leading divergence is the highest of the two possible orders.
    # The leading pole order k is max(k_uv, k_ir).
    print("The leading pole order 'k' is determined by the maximum of UV and IR pole orders.")
    print(f"k_uv = {max_uv_pole_order}")
    print(f"k_ir = {max_ir_pole_order}")

    leading_pole_order = max(max_uv_pole_order, max_ir_pole_order)

    # We print the equation for clarity as requested.
    print(f"Equation for the leading pole order k: max({max_uv_pole_order}, {max_ir_pole_order}) = {leading_pole_order}")

    # The divergent term is proportional to 1/epsilon^k. The question asks for the
    # number of the power, which is the exponent, -k.
    power_of_leading_term = -leading_pole_order

    print("\nThe power of the leading divergent term is:")
    print(power_of_leading_term)
    print("-" * 20)

if __name__ == "__main__":
    solve_phi3_theory_problems()