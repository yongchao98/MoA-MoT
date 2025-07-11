def solve_barbell_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a barbell graph
    with 10 nodes using uniform gossiping.
    """
    # 1. Define the parameters of the barbell graph.
    total_nodes = 10
    # A barbell graph connects two identical cliques with a single bridge edge.
    nodes_per_clique = total_nodes // 2
    num_bridge_edges = 1

    # 2. Calculate the number of edges in one clique.
    # A clique is a complete graph, where the number of edges in a
    # complete graph with 'm' nodes is m * (m - 1) / 2.
    edges_in_clique = nodes_per_clique * (nodes_per_clique - 1) // 2

    # 3. Calculate the total number of edges in the entire graph.
    # Total edges = (edges in clique 1) + (edges in clique 2) + (bridge edge).
    total_edges = edges_in_clique + edges_in_clique + num_bridge_edges

    # 4. In uniform gossiping, every edge has an equal chance of being selected.
    # The probability of sampling the single bottleneck edge is 1 / total_edges.
    probability_of_bottleneck = 1 / total_edges

    # 5. Print the step-by-step calculation.
    print(f"Step 1: Determine the structure of the 10-node barbell graph.")
    print(f"The graph consists of two complete graphs (cliques) of {nodes_per_clique} nodes each, connected by a single bottleneck edge.")
    print("-" * 30)
    
    print(f"Step 2: Calculate the number of edges.")
    print(f"Edges in one K{nodes_per_clique} clique = ({nodes_per_clique} * ({nodes_per_clique} - 1)) / 2 = {edges_in_clique}")
    print(f"Total edges = (edges in clique 1) + (edges in clique 2) + (bottleneck edge)")
    print(f"Total edges = {edges_in_clique} + {edges_in_clique} + {num_bridge_edges} = {total_edges}")
    print("-" * 30)

    print(f"Step 3: Calculate the probability of sampling the bottleneck edge.")
    print("This is 1 divided by the total number of edges.")
    print(f"Final Equation: P(bottleneck) = 1 / {total_edges}")
    print(f"Result: {1} / {total_edges} = {probability_of_bottleneck}")


solve_barbell_probability()