# A script to calculate the probability of sampling the bottleneck edge in a barbell graph.

def solve_gossiping_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a 10-node
    barbell graph under randomized uniform gossiping.
    """
    # Total number of nodes in the graph
    total_nodes = 10

    # A standard 10-node barbell graph consists of two 5-node cliques.
    clique_size = 5

    # The bottleneck edge connects one node from each clique. Let's call them u and v.
    # The degree of such a node is the number of other nodes in its clique (clique_size - 1)
    # plus the one connection to the other clique.
    degree_bottleneck_node = (clique_size - 1) + 1

    # In randomized uniform gossiping, a node is picked uniformly at random.
    # The probability of picking any single node is 1 / total_nodes.
    
    # Then, the chosen node picks a neighbor uniformly at random.
    # The probability of a bottleneck node choosing its counterpart is 1 / its degree.
    
    # The overall probability is the sum of two cases:
    # 1. Node u is picked, and it chooses node v: P(u -> v)
    # 2. Node v is picked, and it chooses node u: P(v -> u)
    # P(total) = P(u -> v) + P(v -> u)
    # P(total) = (1/total_nodes) * (1/degree(u)) + (1/total_nodes) * (1/degree(v))

    prob_u_to_v = (1 / total_nodes) * (1 / degree_bottleneck_node)
    prob_v_to_u = (1 / total_nodes) * (1 / degree_bottleneck_node)

    total_probability = prob_u_to_v + prob_v_to_u

    print("The probability is calculated as: P(pick u) * P(u->v) + P(pick v) * P(v->u)")
    print(f"The final equation with the problem's numbers is:")
    print(f"(1 / {total_nodes}) * (1 / {degree_bottleneck_node}) + (1 / {total_nodes}) * (1 / {degree_bottleneck_node}) = {total_probability}")
    print(f"\nThis simplifies to {prob_u_to_v:.3f} + {prob_v_to_u:.3f} = {total_probability:.3f}.")
    print("As a fraction, the probability is 2/50, which simplifies to 1/25.")


solve_gossiping_probability()
<<<0.04>>>