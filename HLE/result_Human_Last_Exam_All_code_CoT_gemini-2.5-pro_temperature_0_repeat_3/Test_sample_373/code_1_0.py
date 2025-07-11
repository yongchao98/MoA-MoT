def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes using randomized uniform gossiping.
    """
    # Total number of nodes in the barbell graph
    total_nodes = 10

    # A 10-node barbell graph consists of two K5 cliques.
    # The size of each clique is:
    clique_size = total_nodes / 2

    # The bottleneck edge connects a node 'u' from the first clique
    # and a node 'v' from the second clique.

    # Calculate the degree of the node 'u' at one end of the bottleneck.
    # It is connected to (clique_size - 1) nodes in its own clique,
    # plus 1 connection via the bottleneck edge.
    degree_u = (clique_size - 1) + 1

    # By symmetry, the degree of node 'v' is the same.
    degree_v = degree_u

    # The probability of sampling the bottleneck edge (u, v) is the sum of two cases:
    # 1. Node u is picked, and it chooses v.
    # 2. Node v is picked, and it chooses u.
    # Formula: P = (1/N) * (1/deg(u)) + (1/N) * (1/deg(v))

    # Calculate the probability of the first case
    prob_u_to_v = (1 / total_nodes) * (1 / degree_u)

    # Calculate the probability of the second case
    prob_v_to_u = (1 / total_nodes) * (1 / degree_v)

    # The total probability is the sum of the two cases
    total_probability = prob_u_to_v + prob_v_to_u

    # Print the step-by-step calculation
    print("The probability of sampling the bottleneck edge is calculated as follows:")
    print("P(sample edge) = P(select u) * P(u selects v) + P(select v) * P(v selects u)")
    print(f"P(sample edge) = (1/{int(total_nodes)}) * (1/{int(degree_u)}) + (1/{int(total_nodes)}) * (1/{int(degree_v)})")
    print(f"P(sample edge) = {prob_u_to_v} + {prob_v_to_u}")
    print(f"P(sample edge) = {total_probability}")
    print("\nAs a fraction, the probability is 2/50, which simplifies to 1/25.")

solve_gossip_probability()