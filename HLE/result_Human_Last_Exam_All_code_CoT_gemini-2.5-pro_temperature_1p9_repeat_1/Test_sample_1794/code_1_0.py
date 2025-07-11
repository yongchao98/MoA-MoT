def solve_feynman_queries():
    """
    This script calculates the answers to the two physics questions posed by the user.
    """

    # Part 1: Number of distinct planar graphs
    # --------------------------------------------
    # The number of primitive 3-loop 4-point topologies that can be drawn planarly is 2.
    # For the cyclic leg ordering (1, 2, 4, 3), there are 2 possible planar channels
    # (s-channel and u-channel).
    # The total number of distinct planar graphs is the product of these two numbers.

    num_primitive_topologies = 2
    num_planar_channels = 2
    total_distinct_graphs = num_primitive_topologies * num_planar_channels

    print("--- Part 1: Number of Distinct Planar Graphs ---")
    print(f"Number of primitive planar topologies: {num_primitive_topologies}")
    print(f"Number of allowed planar channels for the (1,2,4,3) ordering: {num_planar_channels}")
    print(f"Total number of distinct planar graphs is {num_primitive_topologies} * {num_planar_channels} = {total_distinct_graphs}")
    print("")


    # Part 2: Power of the leading divergent term
    # ----------------------------------------------
    # For massless scattering amplitudes, the leading infrared divergence at L-loops
    # in dimensional regularization (d=4-2epsilon) generally produces a pole of order 1/epsilon^(2L).
    # The number of loops is given as L=3.

    loop_order = 3
    leading_divergence_power = 2 * loop_order

    print("--- Part 2: Power of the Leading Divergence ---")
    print(f"Loop order L: {loop_order}")
    print("The power of the leading 1/epsilon pole for a massless L-loop amplitude is 2*L.")
    print(f"The calculated power for L=3 is 2 * {loop_order} = {leading_divergence_power}")

solve_feynman_queries()
