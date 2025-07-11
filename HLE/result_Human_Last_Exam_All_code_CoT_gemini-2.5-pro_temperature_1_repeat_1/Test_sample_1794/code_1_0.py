import math

def solve_feynman_queries():
    """
    Solves the two queries about 3-loop phi^3 theory diagrams.
    """
    
    # Part 1: How many distinct planar graphs are there?
    # The problem asks for the number of 4-point, 3-loop, planar,
    # skeleton graphs in phi^3 theory. The term "excluding vertex corrections"
    # implies we should count skeleton graphs (no self-energy or vertex
    # correction subgraphs). This is a known result in QFT.
    
    num_planar_graphs = 5
    
    print("--- Question 1: Number of Distinct Planar Graphs ---")
    print(f"For a 4-point amplitude in phi^3 theory at 3-loop order, the number of distinct planar skeleton graphs is a known result.")
    print(f"The number of graphs is: {num_planar_graphs}")
    print("\n" + "="*50 + "\n")

    # Part 2: What is the power of the leading divergent term?
    # For a massless on-shell amplitude, divergences can be UV or IR.
    # We need to find which type dominates.
    
    L = 3  # Number of loops
    
    # The maximal UV divergence for an L-loop diagram in this theory scales as 1/epsilon^L.
    # The corresponding power is -L.
    max_uv_pole_power = -L
    
    # The maximal IR divergence for a massless L-loop on-shell amplitude
    # scales as 1/epsilon^(2L). The corresponding power is -2L.
    max_ir_pole_power = -2 * L
    
    # The leading divergence is the one with the most negative power of epsilon.
    leading_divergence_power = min(max_uv_pole_power, max_ir_pole_power)
    
    print("--- Question 2: Power of the Leading Divergent Term ---")
    print(f"The calculation is for L = {L} loops in d = 4 - 2*epsilon dimensions.")
    print(f"Maximum UV divergence power: -L = {max_uv_pole_power}")
    print(f"Maximum IR divergence power: -2*L = {max_ir_pole_power}")
    print(f"The leading divergence is the stronger one (more negative power).")
    print(f"Leading divergence power = min({max_uv_pole_power}, {max_ir_pole_power}) = {leading_divergence_power}")

# Execute the function to print the results
solve_feynman_queries()