import networkx as nx

def count_biconnected_graphs():
    """
    This function counts the number of 2-vertex-connected simple nonisomorphic
    graphs with 5 vertices.

    It works by iterating through all known nonisomorphic simple connected graphs
    with 5 vertices and checking each one for the 2-vertex-connected property.
    """
    # The networkx `graph_atlas_g` provides a list of all nonisomorphic simple
    # connected graphs with 3 to 7 vertices.
    all_atlas_graphs = nx.graph_atlas_g()

    # We first filter the atlas to get only the graphs with exactly 5 vertices.
    graphs_v5 = [g for g in all_atlas_graphs if g.number_of_vertices() == 5]
    total_connected_graphs = len(graphs_v5)

    # Now, we count how many of these are 2-vertex-connected (biconnected).
    biconnected_count = 0
    for graph in graphs_v5:
        # The `nx.is_biconnected` function checks if a graph is 2-vertex-connected.
        if nx.is_biconnected(graph):
            biconnected_count += 1

    # Calculate the number of graphs that are not biconnected.
    non_biconnected_count = total_connected_graphs - biconnected_count

    # Print the breakdown of the calculation as a simple equation.
    print(f"Total simple nonisomorphic connected graphs with 5 vertices: {total_connected_graphs}")
    print(f"Number of graphs that are NOT 2-vertex-connected: {non_biconnected_count}")
    print("---------------------------------------------------------------------")
    print(f"Number of 2-vertex-connected graphs = {total_connected_graphs} - {non_biconnected_count} = {biconnected_count}")

# Run the function to print the result.
count_biconnected_graphs()
