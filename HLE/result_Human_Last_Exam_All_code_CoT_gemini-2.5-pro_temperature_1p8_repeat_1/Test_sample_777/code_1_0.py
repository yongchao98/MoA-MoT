import networkx as nx
import matplotlib.pyplot as plt

def long_cycle_to_disjoint_cycles_reduction(graph_edges, l):
    """
    Performs a reduction from the Long Cycle problem to the Disjoint Cycles problem.

    Args:
        graph_edges (list of tuples): The edge list for the Long Cycle graph G_L.
        l (int): The required minimum cycle length for Long Cycle.

    Returns:
        A tuple containing the new graph G_D (as an edge list) and the new parameter k.
    """
    print(f"--- Reduction from LongCycle(G_L, l={l}) to DisjointCycles(G_D, k) ---")

    # The parameter k for DisjointCycles is the same as l for LongCycle
    k = l
    print(f"Set new parameter k = l = {k}")

    # The new graph G_D starts with all nodes and edges from G_L
    new_graph_edges = list(graph_edges)
    
    # Find the maximum node index in the original graph to avoid collisions
    max_node_idx = 0
    if graph_edges:
        max_node_idx = max(max(edge) for edge in graph_edges)

    start_node_for_new_cycles = max_node_idx + 1

    # Add k-1 new, disjoint cycles of length k
    print(f"Adding {k-1} new disjoint cycles, each of length {k}.")
    
    current_node = start_node_for_new_cycles
    for i in range(k - 1):
        # Nodes for the i-th new cycle
        # e.g., for cycle 0, nodes are {current, current+1, ..., current+k-1}
        cycle_nodes = list(range(current_node, current_node + k))
        
        # Add edges for this cycle
        for j in range(k):
            u = cycle_nodes[j]
            v = cycle_nodes[(j + 1) % k]
            new_graph_edges.append(tuple(sorted((u, v))))
        
        # Update current_node for the next disjoint cycle
        current_node += k

    print("\n--- Resulting DisjointCycles Instance ---")
    print(f"Parameter k: {k}")
    # Printing a subset of edges for brevity if the list is too long
    edges_to_show = new_graph_edges[:20]
    ellipsis = "..." if len(new_graph_edges) > 20 else ""
    print(f"Graph G_D Edges (showing up to 20): {edges_to_show} {ellipsis}")
    num_nodes_original = max_node_idx + 1
    num_nodes_added = (k-1) * k
    num_nodes_total = num_nodes_original + num_nodes_added
    print(f"Total Vertices: {num_nodes_total} ({num_nodes_original} from G_L + {num_nodes_added} added)")
    print(f"Total Edges: {len(new_graph_edges)}")

    return new_graph_edges, k

if __name__ == '__main__':
    # Example usage:
    # An instance of LongCycle: a 5x5 grid graph and l=10
    # A 5x5 grid has a cycle of length 20, so it has one of length >= 10.
    G = nx.grid_2d_graph(5, 5)
    # Convert nodes from (x,y) tuples to integers for clarity
    node_mapping = {node: i for i, node in enumerate(G.nodes())}
    G_L_edges = [tuple(sorted((node_mapping[u], node_mapping[v]))) for u, v in G.edges()]
    l_param = 10
    
    # Run the reduction
    long_cycle_to_disjoint_cycles_reduction(G_L_edges, l_param)
