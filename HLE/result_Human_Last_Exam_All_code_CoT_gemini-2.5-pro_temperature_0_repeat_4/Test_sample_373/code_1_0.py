def calculate_bottleneck_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a barbell graph
    with 10 nodes during randomized uniform gossiping.
    """
    # Total number of nodes in the graph
    total_nodes = 10

    # A standard 10-node barbell graph consists of two 5-node cliques.
    clique_size = 5

    # The bottleneck edge connects one node from each clique.
    # Let's find the degree of a node connected to the bottleneck edge.
    # It's connected to (clique_size - 1) nodes in its own clique
    # and 1 node (the other end of the bottleneck) in the other clique.
    degree_of_bottleneck_node = (clique_size - 1) + 1

    # The probability of selecting any single node is 1 / total_nodes
    p_select_node = 1 / total_nodes

    # The probability of a selected bottleneck node choosing its neighbor
    # across the bottleneck is 1 / its degree.
    p_choose_bottleneck_neighbor = 1 / degree_of_bottleneck_node

    # The total probability is the sum of two scenarios:
    # 1. Select the first bottleneck node AND it chooses the second.
    # 2. Select the second bottleneck node AND it chooses the first.
    prob_scenario_1 = p_select_node * p_choose_bottleneck_neighbor
    prob_scenario_2 = p_select_node * p_choose_bottleneck_neighbor
    total_probability = prob_scenario_1 + prob_scenario_2

    # Print the equation with the numbers used in the calculation
    print(f"The probability is calculated as: P(select node u) * P(u chooses v) + P(select node v) * P(v chooses u)")
    print(f"({int(1/p_select_node)}/{total_nodes}) * ({int(1/p_choose_bottleneck_neighbor)}/{degree_of_bottleneck_node}) + ({int(1/p_select_node)}/{total_nodes}) * ({int(1/p_choose_bottleneck_neighbor)}/{degree_of_bottleneck_node}) = {total_probability}")

calculate_bottleneck_probability()