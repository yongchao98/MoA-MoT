import fractions

def calculate_bottleneck_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a barbell graph
    with 10 nodes using randomized uniform gossiping.

    The barbell graph is assumed to be two K_5 cliques connected by one edge.
    """

    # 1. Define graph parameters
    total_nodes_N = 10
    # A 10-node barbell graph is typically two K_5 cliques.
    nodes_per_clique_m = 5

    # 2. Calculate the degree of the two nodes connected by the bottleneck.
    # Each bottleneck node is connected to (m-1) nodes in its own clique
    # plus 1 node across the bridge.
    degree_bottleneck_node = (nodes_per_clique_m - 1) + 1

    # 3. Calculate the probability of sampling the bottleneck edge.
    # The probability of sampling an edge (u, v) is the sum of probabilities of
    # two events:
    # a) Picking node u, then its neighbor v.
    # b) Picking node v, then its neighbor u.
    # P(edge) = P(pick u)*P(pick v|u) + P(pick v)*P(pick u|v)
    # P(edge) = (1/N) * (1/degree(u)) + (1/N) * (1/degree(v))
    # Here, degree(u) = degree(v) = degree_bottleneck_node.
    
    # Using fractions for clear representation
    prob_pick_node = fractions.Fraction(1, total_nodes_N)
    prob_pick_neighbor = fractions.Fraction(1, degree_bottleneck_node)
    
    total_prob = prob_pick_node * prob_pick_neighbor + prob_pick_node * prob_pick_neighbor

    # 4. Print the explanation and the final equation.
    print("A barbell graph with 10 nodes consists of two 5-node cliques (K_5) connected by a single bottleneck edge.")
    print("In randomized uniform gossiping, a node is picked uniformly at random, then one of its neighbors is picked uniformly at random.")
    print("\nLet's calculate the probability of sampling the bottleneck edge.")
    print(f"Total number of nodes (N): {total_nodes_N}")
    print(f"Degree of each of the two nodes on the bottleneck edge (d): (5-1) + 1 = {degree_bottleneck_node}")
    
    print("\nThe probability equation is: P = (1/N) * (1/d) + (1/N) * (1/d)")
    
    # Printing the final equation with the numbers plugged in
    print("\nFinal Equation:")
    print(f"P = (1 / {total_nodes_N}) * (1 / {degree_bottleneck_node}) + (1 / {total_nodes_N}) * (1 / {degree_bottleneck_node})")
    print(f"P = {prob_pick_node * prob_pick_neighbor} + {prob_pick_node * prob_pick_neighbor}")
    print(f"P = {total_prob}")
    print(f"\nThe final probability is {float(total_prob):.2f}.")

calculate_bottleneck_probability()
<<<0.04>>>