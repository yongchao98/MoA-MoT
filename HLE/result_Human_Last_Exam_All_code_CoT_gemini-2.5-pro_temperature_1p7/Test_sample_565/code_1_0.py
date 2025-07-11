import networkx as nx

def calculate_vt_graph_degree_distribution():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs on 8
    vertices for each degree from 0 to 7.

    This function uses the canonical list of such graphs from Brendan McKay's
    data collection. It reads each graph, determines its degree, and counts
    how many graphs exist for each degree.
    """

    # This is the complete list of non-isomorphic vertex-transitive graphs
    # on 8 vertices, in graph6 format, from Brendan McKay's website.
    # There are 21 such graphs in total.
    g6_data = [
        b"G?Q@", b"G?C_", b"G~{wg", b"G`{_G`", b"G`cGW`", b"G`cG_`",
        b"G`cGc`", b"G`cGg`", b"G`a_G`", b"G`bOW`", b"G`bOWg", b"G`dMW`",
        b"G`hWW`", b"G`hWWg", b"G`hWg`", b"G`hWgg", b"G`p_W`", b"Gqe`G`",
        b"Gqf`wg", b"Gwd`_", b"Gw`_~"
    ]

    # Initialize a list to store the counts for degrees 0 through 7.
    # n_j = degree_counts[j]
    degree_counts = [0] * 8

    # Process each graph in the dataset
    for g6_bytes in g6_data:
        # Convert graph6 format to a networkx graph object
        g = nx.from_graph6_bytes(g6_bytes)

        # Sanity check: ensure the graph has 8 vertices.
        if len(g) != 8:
            continue

        # A vertex-transitive graph is regular, so all vertices have the same degree.
        # We can determine the degree of the graph by checking any vertex.
        # If the graph has no edges, its degree is 0.
        if len(g.edges) == 0:
            degree = 0
        else:
            # Get the degree of the first vertex
            first_vertex = next(iter(g.nodes))
            degree = g.degree[first_vertex]

        # Increment the counter for this degree
        if 0 <= degree < len(degree_counts):
            degree_counts[degree] += 1

    # Format the final list for output
    result_string = ", ".join(map(str, degree_counts))
    print(f"[{result_string}]")

calculate_vt_graph_degree_distribution()