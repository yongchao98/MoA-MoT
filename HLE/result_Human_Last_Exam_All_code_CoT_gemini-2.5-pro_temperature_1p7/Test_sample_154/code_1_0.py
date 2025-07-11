import networkx as nx

def count_2_vertex_connected_graphs():
    """
    Calculates the number of 2-vertex-connected simple nonisomorphic graphs
    with a specified number of vertices.
    """
    num_vertices = 5
    count = 0

    # Generate all simple nonisomorphic graphs with 5 vertices.
    # The networkx generator handles the "simple" and "nonisomorphic" constraints.
    non_iso_graphs = nx.nonisomorphic_graphs(num_vertices)

    # Iterate through each unique graph
    for G in non_iso_graphs:
        # A graph is 2-vertex-connected if its node connectivity is at least 2.
        # This means the graph cannot be disconnected by removing just one vertex.
        if nx.node_connectivity(G) >= 2:
            count += 1
            
    # The problem asks to "output each number in the final equation".
    # We will print a descriptive sentence including the input number (5) and the result.
    print(f"Number of vertices: {num_vertices}")
    print(f"Number of 2-vertex-connected simple nonisomorphic graphs found: {count}")

count_2_vertex_connected_graphs()