import math

def solve_barbell_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a 10-node barbell graph.
    """
    # Total number of nodes in the barbell graph
    total_nodes = 10

    # A barbell graph with 10 nodes consists of two equal-sized cliques
    # connected by a single edge (the bottleneck).
    nodes_per_clique = total_nodes // 2

    # The number of edges in a single complete graph (clique) with 'n' nodes
    # is given by the formula n * (n - 1) / 2.
    edges_in_one_clique = nodes_per_clique * (nodes_per_clique - 1) // 2

    # The total number of edges is the sum of edges in both cliques plus the
    # one bottleneck edge.
    num_bottleneck_edges = 1
    total_edges = (2 * edges_in_one_clique) + num_bottleneck_edges

    # In uniform randomized gossiping, each edge is chosen with equal probability.
    # The probability of sampling the bottleneck edge is 1 / (total number of edges).
    probability = 1 / total_edges

    print(f"A barbell graph with {total_nodes} nodes is composed of two K{nodes_per_clique} cliques.")
    print(f"Number of edges in one K{nodes_per_clique} clique = ({nodes_per_clique} * ({nodes_per_clique} - 1)) / 2 = {edges_in_one_clique}")
    print(f"Total number of edges in the two cliques = 2 * {edges_in_one_clique} = {2 * edges_in_one_clique}")
    print(f"Number of bottleneck edges connecting the cliques = {num_bottleneck_edges}")
    print(f"Total edges in the graph = {2 * edges_in_one_clique} + {num_bottleneck_edges} = {total_edges}")
    print("\nIn randomized uniform gossiping, the probability of sampling a specific edge is 1 / (Total Edges).")
    print(f"The probability of sampling the bottleneck edge is: 1 / {total_edges}")
    print(f"Calculated probability: {probability:.4f}")

solve_barbell_probability()