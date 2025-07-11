import networkx as nx

def solve():
    """
    Calculates the number of isomorphism classes of vertex-transitive graphs
    with 8 vertices for each possible vertex degree j, from j=0 to j=7.
    """
    N = 8
    # Initialize a list to store the counts for each degree j=0..7
    n_counts = [0] * N

    # The graph_atlas_g() function provides all non-isomorphic graphs.
    # We will iterate through them and filter for those with 8 vertices.
    all_graphs = nx.graph_atlas_g()
    
    graphs_n8 = []
    for G in all_graphs:
        if G.number_of_nodes() == N:
            graphs_n8.append(G)

    # For each non-isomorphic graph with 8 vertices, check its properties.
    for G in graphs_n8:
        # A vertex-transitive graph must be regular (all vertices have the same degree).
        # First, check for regularity.
        degrees = [d for v, d in G.degree()]
        if not degrees:  # Handle the case of a graph with 0 vertices if it were possible
            continue
            
        first_degree = degrees[0]
        is_regular = all(d == first_degree for d in degrees)

        if is_regular:
            # If the graph is regular, check for vertex-transitivity.
            if nx.is_vertex_transitive(G):
                # If it is, increment the counter for its degree.
                degree_j = first_degree
                if 0 <= degree_j < N:
                    n_counts[degree_j] += 1
    
    # Print the final list of counts in the required format.
    print(n_counts)

solve()