def calculate_bottleneck_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes during randomized uniform gossiping.
    """
    # Define the parameters of the barbell graph
    total_nodes = 10
    num_cliques = 2
    nodes_per_clique = total_nodes // num_cliques

    # Calculate the degree of one of the two nodes connected by the bottleneck edge.
    # Each bottleneck node is connected to (nodes_per_clique - 1) nodes
    # in its own clique, plus the other bottleneck node.
    degree_of_bottleneck_node = (nodes_per_clique - 1) + 1

    # In uniform gossiping, any node is chosen with probability 1 / total_nodes.
    prob_selecting_a_node = f"1/{total_nodes}"

    # Given a bottleneck node is chosen, the probability of it selecting the
    # other bottleneck node (its neighbor across the bridge) is 1 / its degree.
    prob_node_chooses_bottleneck_neighbor = f"1/{degree_of_bottleneck_node}"

    # The total probability is the sum of two scenarios:
    # 1. First bottleneck node is chosen AND it selects the second.
    # 2. Second bottleneck node is chosen AND it selects the first.
    # P = P(select node u) * P(u selects v) + P(select node v) * P(v selects u)
    prob_scenario_1 = (1 / total_nodes) * (1 / degree_of_bottleneck_node)
    prob_scenario_2 = (1 / total_nodes) * (1 / degree_of_bottleneck_node)
    total_probability = prob_scenario_1 + prob_scenario_2

    print("The probability of sampling the bottleneck edge is calculated as follows:")
    print("P = P(select first bottleneck node) * P(it chooses the bottleneck edge) + P(select second bottleneck node) * P(it chooses the bottleneck edge)")
    print("\nThe final equation is:")
    print(f"P = ({prob_selecting_a_node}) * ({prob_node_chooses_bottleneck_neighbor}) + ({prob_selecting_a_node}) * ({prob_node_chooses_bottleneck_neighbor})")
    print(f"P = {prob_scenario_1} + {prob_scenario_2}")
    print(f"P = {total_probability}")


calculate_bottleneck_gossip_probability()