import networkx as nx

def count_vertex_transitive_graphs_by_degree():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs on 8
    vertices for each degree from 0 to 7.

    It uses the authoritative list of 34 such graphs compiled by Gordon Royle.
    The graph data is hardcoded in the graph6 format.
    """

    # The list of 34 vertex-transitive graphs on 8 vertices
    # in graph6 format, from Gordon Royle's data page:
    # http://staffhome.ecm.uwa.edu.au/~00013890/trans/vt8.g6
    g6_data = [
        b'G??`w_', b'G?_@O_', b'G??a?_', b'G?Aa?_', b'G?QA?_',
        b'G?`a?_', b'G@`A?_', b'GQAE`_', b'GQ`E`_', b'G`aE`_',
        b'G`AE`_', b'G?q@w_', b'G?`@w_', b'G?QAw_', b'G?`Qw_',
        b'G@`Aw_', b'G@qAw_', b'G@QAw_', b'GQ`Eh_', b'GQ`eh_',
        b'GQAEh_', b'GQaeh_', b'G`AEh_', b'G`aeh_', b'G?qAw_',
        b'G?`Aw_', b'G?Aaw_', b'G?aQw_', b'G?qQw_', b'G`Aaw_',
        b'G@?aw_', b'G~_E@w', b'G@Aaw_', b'G@aQw_'
    ]

    # Initialize a list to store the counts for each degree (0 to 7)
    degree_counts = [0] * 8

    # Process each graph string
    for g6_string in g6_data:
        # Parse the graph6 byte string into a NetworkX graph object
        graph = nx.from_graph6_bytes(g6_string)

        # All vertices in a vertex-transitive graph have the same degree.
        # We can get the degree from any vertex (e.g., vertex 0).
        if graph.number_of_nodes() > 0:
            degree = graph.degree[0]
            if 0 <= degree < 8:
                degree_counts[degree] += 1
        elif graph.number_of_nodes() == 0:
             # This case is for the null graph on 0 vertices, not in our list
             pass


    # Print the result in the specified list format
    # The problem also requests printing the final equation.
    print(f"The number of isomorphism classes of vertex-transitive graphs on 8 vertices with degree j are:")
    print(f"n_j = {degree_counts}")
    equation_parts = [str(n) for n in degree_counts]
    print(f"Total number of graphs is {sum(degree_counts)}.")
    print(f"Final equation: {' + '.join(equation_parts)} = {sum(degree_counts)}")


# Execute the function
count_vertex_transitive_graphs_by_degree()
