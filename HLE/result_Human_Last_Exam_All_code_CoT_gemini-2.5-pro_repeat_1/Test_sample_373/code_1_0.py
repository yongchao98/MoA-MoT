def solve_gossiping_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes using randomized uniform gossiping.
    """
    # Total number of nodes in the barbell graph
    total_nodes = 10

    # The barbell graph has two cliques. The number of nodes in each clique is:
    nodes_per_clique = total_nodes / 2

    # Let the two nodes connected by the bottleneck edge be 'u' and 'v'.
    # A node's degree is its number of connections. In a clique of size k, a node
    # connects to k-1 other nodes. The bottleneck node also has 1 connection
    # to the other clique.
    # Degree = (nodes_per_clique - 1) + 1
    degree_bottleneck_node = nodes_per_clique

    # In uniform gossiping, a node is first chosen uniformly at random.
    # The probability of picking any single node is 1 / total_nodes.
    prob_pick_node = 1 / total_nodes

    # If we pick a node, the probability of then picking a specific neighbor
    # is 1 / degree of that node.
    prob_pick_neighbor = 1 / degree_bottleneck_node

    # The bottleneck edge can be sampled in two ways:
    # 1. Pick node 'u', then its neighbor 'v' across the bottleneck.
    prob_u_then_v = prob_pick_node * prob_pick_neighbor
    # 2. Pick node 'v', then its neighbor 'u' across the bottleneck.
    prob_v_then_u = prob_pick_node * prob_pick_neighbor

    # The total probability is the sum of these two mutually exclusive events.
    total_probability = prob_u_then_v + prob_v_then_u

    print("Step 1: Define the graph properties.")
    print(f"Total nodes (N): {total_nodes}")
    print(f"Degree of each node on the bottleneck edge (deg): {int(degree_bottleneck_node)}")
    print("-" * 30)

    print("Step 2: Define the probability formula for sampling edge (u, v).")
    print("P(sample edge) = P(pick u) * P(pick v | u) + P(pick v) * P(pick u | v)")
    print("P(sample edge) = (1/N) * (1/deg(u)) + (1/N) * (1/deg(v))")
    print("-" * 30)

    print("Step 3: Substitute the values into the formula.")
    print(f"P(sample edge) = (1/{total_nodes}) * (1/{int(degree_bottleneck_node)}) + (1/{total_nodes}) * (1/{int(degree_bottleneck_node)})")
    print(f"P(sample edge) = {prob_u_then_v} + {prob_v_then_u}")
    print("-" * 30)
    
    print("Final Answer:")
    print(f"The total probability of sampling the bottleneck edge is: {total_probability}")

solve_gossiping_probability()