def solve_barbell_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a barbell graph.

    A barbell graph with 10 nodes is defined as two complete 5-node graphs (K5)
    connected by a single bridge edge (the bottleneck).

    Randomized uniform gossiping involves selecting one edge uniformly at random
    from all edges in the graph.
    """
    # Parameters based on the problem description
    total_nodes = 10
    num_cliques = 2
    num_bottleneck_edges = 1

    # 1. Calculate the number of nodes per clique
    nodes_per_clique = total_nodes // num_cliques

    # 2. Calculate the number of edges in one K5 clique
    # The formula for edges in a complete graph with n nodes is n*(n-1)/2
    edges_per_clique = nodes_per_clique * (nodes_per_clique - 1) // 2

    # 3. The final probability is the ratio of bottleneck edges to the total number of edges.
    # The total number of edges is the sum of edges from both cliques plus the bottleneck edge.
    # We will show the components of this sum in the final output.
    numerator = num_bottleneck_edges
    edges_clique1 = edges_per_clique
    edges_clique2 = edges_per_clique
    
    # To satisfy the request to "output each number in the final equation",
    # we will print the probability as a fraction showing how the denominator is composed.
    print(f"{numerator} / ({edges_clique1} + {edges_clique2} + {num_bottleneck_edges})")

solve_barbell_probability()