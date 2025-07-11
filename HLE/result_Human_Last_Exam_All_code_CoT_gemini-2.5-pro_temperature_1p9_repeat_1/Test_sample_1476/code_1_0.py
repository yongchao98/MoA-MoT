import networkx as nx
import numpy as np

def illustrate_graph_signal_inference():
    """
    This function illustrates the inference based on the problem description.
    1. Create a sample graph.
    2. Define a vertex signal x0 that must be constant on connected components
       to satisfy the problem's premises.
    3. Compute the edge signal x1 from x0.
    4. Verify that the conditions are met.
    5. Compute and print the total variation, showing it's zero.
    """
    # 1. Create a sample graph (e.g., a lollipop graph, which has a cycle)
    G = nx.lollipop_graph(4, 2) # A K4 clique attached to a path of length 2
    nodes = list(G.nodes())
    edges = list(G.edges())
    num_nodes = len(nodes)
    num_edges = len(edges)
    
    print(f"Graph has {num_nodes} nodes and {num_edges} edges.")
    print("Nodes:", nodes)
    print("Edges:", edges)
    
    # The premises force x0 to be constant on connected components.
    # Since the graph is connected, x0 must be a constant vector.
    constant_value = 10.0
    x0 = np.full(num_nodes, constant_value)
    print("\nVertex signal x0:", x0)

    # 3. Compute the edge signal x1 using the given formula
    # x1_e = |x0_u - x0_v|
    x1 = np.zeros(num_edges)
    edge_map = {edge: i for i, edge in enumerate(edges)}
    for u, v in edges:
        idx = edge_map[(u, v)]
        x1[idx] = np.abs(x0[u] - x0[v])

    print("Edge signal x1:", x1)
    
    # 4. Verification of premises (conceptual)
    # Premise 1: Sum over cycles is 0. Since x1 is all zeros, this holds.
    # Premise 2: Divergence is 0 (B1 * x1 = 0). Since x1 is all zeros, this holds.
    # Premise 3: The relation between x0 and x1 holds by construction.
    
    # 5. Compute and print the Total Variation of x0
    print("\nCalculating Total Variation (TV) of the vertex signal x0:")
    tv_sum_str_parts = []
    tv_val_str_parts = []
    tv_terms = []

    for u, v in edges:
        val = np.abs(x0[u] - x0[v])
        tv_terms.append(val)
        tv_sum_str_parts.append(f"|x0[{u}]-x0[{v}]|")
        tv_val_str_parts.append(f"|{x0[u]}-{x0[v]}|")

    tv = np.sum(tv_terms)

    # Outputting each number in the final equation as requested.
    print(f"TV = {' + '.join(tv_sum_str_parts)}")
    print(f"   = {' + '.join(tv_val_str_parts)}")
    print(f"   = {' + '.join([str(term) for term in tv_terms])}")
    print(f"   = {tv}")

    print("\nConclusion: The premises together imply that the total variation must be 0.")

if __name__ == '__main__':
    illustrate_graph_signal_inference()