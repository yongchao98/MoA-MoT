import math

def solve_barbell_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes.
    """
    total_nodes = 10
    num_cliques = 2
    num_bottleneck_edges = 1

    # For a symmetric barbell graph, the nodes are split evenly between the cliques.
    clique_size = total_nodes // num_cliques

    # The number of edges in a complete graph (clique) with 'm' nodes is m * (m-1) / 2.
    edges_per_clique = (clique_size * (clique_size - 1)) / 2
    edges_per_clique = int(edges_per_clique)

    # Total edges = (edges in clique 1) + (edges in clique 2) + (bottleneck edge).
    total_edges = (num_cliques * edges_per_clique) + num_bottleneck_edges

    # The probability is the ratio of bottleneck edges to total edges.
    probability = num_bottleneck_edges / total_edges

    print(f"A barbell graph with {total_nodes} nodes consists of {num_cliques} cliques of size {clique_size}.")
    print(f"Number of edges in one K{clique_size} clique = ({clique_size} * ({clique_size} - 1)) / 2 = {edges_per_clique}")
    print(f"Total number of edges = (2 * {edges_per_clique}) + {num_bottleneck_edges} = {total_edges}")
    print("\nThe probability of sampling the single bottleneck edge is the ratio of bottleneck edges to total edges.")
    print(f"\nFinal Equation:")
    print(f"P(bottleneck) = {num_bottleneck_edges} / {total_edges}")
    print(f"\nResult: {probability:.4f}")

solve_barbell_probability()
