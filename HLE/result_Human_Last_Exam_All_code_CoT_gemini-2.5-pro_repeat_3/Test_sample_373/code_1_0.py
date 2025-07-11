def solve_barbell_gossip_prob():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes under uniform gossiping.
    """
    # Total number of nodes in the barbell graph
    n = 10

    # A barbell graph with n nodes is formed by two complete graphs (cliques)
    # connected by a single edge (the bottleneck).
    # Each clique has m = n / 2 nodes.
    m = n / 2

    # In a complete graph K_m, each node has a degree of m-1.
    # The two nodes connected by the bottleneck edge have one additional connection.
    # Therefore, their degree is (m - 1) + 1 = m.
    deg_bottleneck_node = m

    print("Step-by-step calculation:")
    print(f"1. The graph is a barbell graph with a total of n = {n} nodes.")
    print(f"2. It consists of two {int(m)}-node complete graphs (cliques) connected by one bottleneck edge.")
    print(f"3. The two nodes connected by the bottleneck edge each have a degree of ({int(m)} - 1) + 1 = {int(deg_bottleneck_node)}.")
    print("\nRandomized uniform gossiping samples an edge by:")
    print("  a) Picking a node 'u' uniformly at random (probability 1/n).")
    print("  b) Picking one of its neighbors 'v' uniformly at random (probability 1/deg(u)).")
    
    print("\nTo sample the bottleneck edge, we must pick one of the two nodes it connects and then its partner.")
    print("Let the bottleneck edge connect nodes 'u_b' and 'v_b'.")
    print(f"We have n = {n}, deg(u_b) = {int(deg_bottleneck_node)}, and deg(v_b) = {int(deg_bottleneck_node)}.")

    # The probability is the sum of two scenarios:
    # 1. Pick the first bottleneck node, then its partner across the bridge.
    # 2. Pick the second bottleneck node, then its partner across the bridge.
    prob_one_way = (1 / n) * (1 / deg_bottleneck_node)
    total_prob = prob_one_way * 2

    # Final output showing the equation
    print("\nThe probability is the sum of sampling the edge from either direction:")
    print(f"P(bottleneck) = P(pick u_b, then v_b) + P(pick v_b, then u_b)")
    print(f"P(bottleneck) = (1/n * 1/deg(u_b)) + (1/n * 1/deg(v_b))")
    print(f"P(bottleneck) = (1/{n} * 1/{int(deg_bottleneck_node)}) + (1/{n} * 1/{int(deg_bottleneck_node)})")
    print(f"P(bottleneck) = {1/n * 1/deg_bottleneck_node} + {1/n * 1/deg_bottleneck_node}")
    print(f"P(bottleneck) = {total_prob}")

solve_barbell_gossip_prob()