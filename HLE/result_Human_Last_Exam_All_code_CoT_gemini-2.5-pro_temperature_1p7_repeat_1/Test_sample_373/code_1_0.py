def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    10-node barbell graph during randomized uniform gossiping.
    """
    # Total number of nodes in the graph
    total_nodes = 10

    # The graph consists of two 5-node complete graphs (K5).
    # We consider the two nodes that form the bottleneck edge. Let's call them u and v.
    
    # Calculate the degree of the first bottleneck node 'u'.
    # It is connected to 4 other nodes in its own K5 clique and 1 node 'v' in the other clique.
    nodes_in_clique = 5
    degree_u = (nodes_in_clique - 1) + 1
    
    # By symmetry, the degree of the second bottleneck node 'v' is the same.
    degree_v = (nodes_in_clique - 1) + 1

    # The probability is the sum of two scenarios:
    # 1. Node 'u' is picked and gossips to 'v': P(pick u) * P(u gossips to v)
    #    = (1 / total_nodes) * (1 / degree_u)
    # 2. Node 'v' is picked and gossips to 'u': P(pick v) * P(v gossips to u)
    #    = (1 / total_nodes) * (1 / degree_v)
    
    # Since these are mutually exclusive events, we sum their probabilities.
    # Total P = (1 / (total_nodes * degree_u)) + (1 / (total_nodes * degree_v))
    
    prob_u_to_v_denominator = total_nodes * degree_u
    prob_v_to_u_denominator = total_nodes * degree_v
    
    # P = 1/50 + 1/50 = 2/50
    total_prob_numerator = 2
    total_prob_denominator = prob_u_to_v_denominator
    
    final_decimal_prob = total_prob_numerator / total_prob_denominator

    # Output the final equation with all numbers
    print(
        f"The probability of sampling the bottleneck edge is calculated by summing the probabilities "
        f"of each bottleneck node being chosen and gossiping to the other."
    )
    print(
        f"The calculation is: (1/{total_nodes}) * (1/{degree_u}) + (1/{total_nodes}) * (1/{degree_v}) "
        f"= 1/{prob_u_to_v_denominator} + 1/{prob_v_to_u_denominator} "
        f"= {total_prob_numerator}/{total_prob_denominator} = {final_decimal_prob}"
    )

solve_gossip_probability()