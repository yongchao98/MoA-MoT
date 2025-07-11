def calculate_bottleneck_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a barbell graph
    with 10 nodes using randomized uniform gossiping.
    """
    # Total number of nodes in the barbell graph
    total_nodes = 10

    # A 10-node barbell graph consists of two K5 cliques.
    # Number of nodes in each clique
    nodes_per_clique = 5

    # The bottleneck edge connects one node from each clique.
    # The degree of such a node is its neighbors within the clique (nodes_per_clique - 1)
    # plus the one connection via the bottleneck.
    degree_bottleneck_node = (nodes_per_clique - 1) + 1

    # In uniform gossiping, the probability of sampling an edge (u, v) is:
    # P(pick u) * P(pick v | u) + P(pick v) * P(pick u | v)
    # P(pick node) = 1 / total_nodes
    # P(pick neighbor | node) = 1 / degree(node)

    prob_pick_u = f"1/{total_nodes}"
    prob_pick_v_from_u = f"1/{degree_bottleneck_node}"
    
    prob_pick_v = f"1/{total_nodes}"
    prob_pick_u_from_v = f"1/{degree_bottleneck_node}"

    # The final probability calculation
    final_prob = (1 / total_nodes) * (1 / degree_bottleneck_node) + \
                 (1 / total_nodes) * (1 / degree_bottleneck_node)
                 
    term1_denominator = total_nodes * degree_bottleneck_node
    term2_denominator = total_nodes * degree_bottleneck_node
    
    final_numerator = 2
    final_denominator = total_nodes * degree_bottleneck_node
    
    simplified_denominator = final_denominator // final_numerator
    
    print("The probability of sampling the bottleneck edge is calculated as follows:")
    print("P(sampling) = P(picking first node) * P(picking second node) + P(picking second node) * P(picking first node)")
    print(f"P(sampling) = ({prob_pick_u}) * ({prob_pick_v_from_u}) + ({prob_pick_v}) * ({prob_pick_u_from_v})")
    print(f"P(sampling) = 1/{term1_denominator} + 1/{term2_denominator}")
    print(f"P(sampling) = {final_numerator}/{final_denominator}")
    print(f"P(sampling) = 1/{simplified_denominator}")
    print(f"\nThe final probability is: {final_prob}")

calculate_bottleneck_probability()