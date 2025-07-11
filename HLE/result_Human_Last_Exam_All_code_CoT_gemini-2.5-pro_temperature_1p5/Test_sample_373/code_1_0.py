def solve_barbell_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a 10-node barbell graph.
    """
    
    # A 10-node barbell graph consists of two cliques.
    # The total number of nodes is 10, so each clique will have 10 / 2 = 5 nodes.
    nodes_per_clique = 5
    
    # The number of edges in a complete graph (clique) with 'n' nodes is n * (n - 1) / 2.
    edges_in_one_clique = nodes_per_clique * (nodes_per_clique - 1) // 2
    
    # The total number of edges is the sum of edges in the two cliques plus the single bottleneck edge connecting them.
    num_cliques = 2
    bottleneck_edges = 1
    total_edges = (num_cliques * edges_in_one_clique) + bottleneck_edges
    
    # The probability of sampling the bottleneck edge is 1 / total_edges.
    numerator = 1
    denominator = total_edges
    probability = numerator / denominator

    print(f"The 10-node barbell graph has two {nodes_per_clique}-node cliques connected by {bottleneck_edges} bridge edge.")
    print(f"The number of edges in one {nodes_per_clique}-node clique is: {edges_in_one_clique}")
    print(f"The total number of edges in the graph is: (2 * {edges_in_one_clique}) + {bottleneck_edges} = {total_edges}")
    print("\nIn randomized uniform gossiping, an edge is chosen uniformly at random.")
    print("The probability of sampling the bottleneck edge is the ratio of bottleneck edges to total edges.")
    print("The final equation is:")
    print(f"Probability = {numerator} / {denominator}")
    print(f"\nThe probability is approximately: {probability:.5f}")

solve_barbell_probability()