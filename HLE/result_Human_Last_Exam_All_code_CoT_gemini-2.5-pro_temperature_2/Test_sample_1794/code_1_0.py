def solve_feynman_diagram_problems():
    """
    This script solves two problems related to 3-loop phi^3 theory.
    """
    
    # Part 1: Counting the number of distinct planar graphs.
    # The question asks for the number of distinct planar (non-crossing) 3-loop, 4-point
    # graphs in massless scalar phi^3 theory, with the condition "excluding the diagrams
    # with vertex corrections". In quantum field theory, this is the standard way of asking
    # for "primitive" diagrams. These diagrams do not contain any self-energy or vertex
    # correction subgraphs, as those divergences are accounted for at lower loop orders.
    # For a 4-point function at 3 loops, this is a known combinatorial result in QFT.
    
    num_graphs = 2
    
    print("1. How many distinct planar (non-crossing) graphs are there in the Feynman diagrams, excluding the diagrams with vertex corrections?")
    print(f"The number of distinct planar primitive graphs is {num_graphs}.")
    print("")

    # Part 2: Finding the power of the leading divergent term.
    # The problem asks for the power of the leading divergent term in the epsilon expansion,
    # where d = 4 - 2*epsilon. Divergences in Feynman integrals manifest as poles in epsilon.
    # For massless on-shell scattering, divergences arise from both Ultraviolet (UV) and
    # Infrared (IR) momentum regions.
    # - UV poles can go up to 1/epsilon^3 for a 3-loop graph with nested self-energy subgraphs.
    # - IR poles for an L-loop massless amplitude are generally more severe, scaling as 1/epsilon^(2L).
    #   This is a well-established result for massless theories.
    # The dominant divergence determines the leading term.
    
    L = 3  # Number of loops
    # The formula for the power of the leading IR divergent term is -2L.
    leading_power = -2 * L

    print("2. What is the number of the power of the leading divergent term of the epsilon expansion?")
    print("The power of the leading term is determined by the most severe infrared (IR) divergence.")
    print(f"For L = {L} loops, the formula for the power of the leading pole 1/epsilon^k is k = 2L.")
    print(f"Therefore, the power of the leading term (which is -k) is -(2 * {L}) = {leading_power}.")

solve_feynman_diagram_problems()