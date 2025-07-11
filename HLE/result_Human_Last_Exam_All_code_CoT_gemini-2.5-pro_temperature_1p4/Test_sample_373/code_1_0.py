def solve_barbell_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes during uniform gossiping.
    """
    # Total number of nodes in the graph
    n = 10

    # A standard barbell graph with n nodes has two cliques of size m=n/2
    m = n / 2

    # The nodes connected by the bottleneck edge (let's call them b1 and b2)
    # are each connected to (m-1) nodes in their own clique plus 1 node
    # across the bridge.
    # So, their degree is (m-1) + 1 = m.
    degree_bridge_node = m

    # The probability of sampling the bottleneck edge (b1, b2) is the sum
    # of two probabilities:
    # 1. Selecting b1 (prob 1/n) AND b1 choosing b2 (prob 1/degree(b1))
    # 2. Selecting b2 (prob 1/n) AND b2 choosing b1 (prob 1/degree(b2))

    # Calculate the probability of one of these symmetric events
    prob_one_direction = (1/n) * (1/degree_bridge_node)

    # The total probability is the sum of the two events
    total_probability = prob_one_direction + prob_one_direction

    print("To find the probability of sampling the bottleneck edge, we follow these steps:")
    print(f"1. A barbell graph with {n} nodes has two cliques of size {int(m)}.")
    print(f"2. The two nodes on the bottleneck edge each have a degree of {int(degree_bridge_node)} ( {int(m-1)} in-clique neighbors + 1 bottleneck neighbor).")
    print(f"3. In uniform gossiping, any node is chosen with probability 1/{n}.")
    print("4. The probability of sampling the bottleneck is the sum of two mutually exclusive events: node 'b1' choosing 'b2', or node 'b2' choosing 'b1'.")
    print("\nThe final equation is:")
    print(f"P(bottleneck) = P(pick b1) * P(b1 picks b2) + P(pick b2) * P(b2 picks b1)")
    print(f"P(bottleneck) = (1/{int(n)}) * (1/{int(degree_bridge_node)}) + (1/{int(n)}) * (1/{int(degree_bridge_node)})")
    print(f"P(bottleneck) = {prob_one_direction} + {prob_one_direction}")
    print(f"P(bottleneck) = {total_probability}")

solve_barbell_gossip_probability()