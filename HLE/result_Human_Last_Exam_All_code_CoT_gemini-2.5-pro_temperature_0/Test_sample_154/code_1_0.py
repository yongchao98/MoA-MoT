import networkx as nx

def solve_graph_counting():
    """
    This function calculates the number of 2-vertex-connected simple
    nonisomorphic graphs with 5 vertices.
    """
    # The networkx graph atlas contains all nonisomorphic graphs up to 7 vertices.
    # We load the entire atlas.
    all_graphs = nx.graph_atlas_g()

    # Filter the list to get only graphs with 5 vertices.
    graphs_with_5_vertices = [g for g in all_graphs if g.number_of_nodes() == 5]

    # Initialize a counter for the graphs that meet the connectivity criteria.
    biconnected_graph_count = 0

    # Iterate through each 5-vertex graph and check for 2-vertex-connectivity.
    for g in graphs_with_5_vertices:
        # A graph is k-vertex-connected if it has more than k vertices and
        # remains connected when fewer than k vertices are removed.
        # nx.is_k_vertex_connected(G, 2) returns True if the graph is 2-vertex-connected.
        if nx.is_k_vertex_connected(g, 2):
            biconnected_graph_count += 1

    # Print the final result.
    print(f"Total number of simple nonisomorphic graphs with 5 vertices: {len(graphs_with_5_vertices)}")
    print(f"Number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices: {biconnected_graph_count}")

if __name__ == "__main__":
    solve_graph_counting()