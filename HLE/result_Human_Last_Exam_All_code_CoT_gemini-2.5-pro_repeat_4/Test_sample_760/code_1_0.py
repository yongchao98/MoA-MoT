import networkx as nx
import numpy as np

def check_graph_fixed_point(G):
    """
    Checks if a graph G is a fixed point of the transformation T.
    A graph is a fixed point if T(G) = G.
    T(G) has an edge (x,y) iff G has 1 or 2 paths of length 2 between x and y.
    """
    n = G.number_of_vertices()
    
    # The case n=0 is not connected.
    if n == 0:
        return False

    # For n=1, the graph is a single vertex. It is connected.
    # The definition of T applies to "distinct vertices x and y".
    # As there are no such pairs, T(G) has no edges, just like G. So T(G) = G.
    if n == 1:
        return True
    
    # For n > 1, we check the condition using adjacency matrices.
    # The number of length-2 paths between vertices is given by A^2.
    try:
        # Use a sorted nodelist to ensure consistent matrix representation.
        nodelist = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodelist)
    except Exception:
        # In case of any issue with graph conversion.
        return False

    A_squared = np.dot(A, A)

    # We need to check if the graph G satisfies the condition for being a fixed point.
    # The condition is: an edge (i,j) exists in G if and only if the number of
    # common neighbors is 1 or 2.
    for i in range(n):
        for j in range(i + 1, n):
            has_edge = A[i, j] == 1
            num_common_neighbors = A_squared[i, j]
            
            # If there is an edge, the number of common neighbors must be 1 or 2.
            if has_edge and num_common_neighbors not in [1, 2]:
                return False
            # If there is no edge, the number of common neighbors must NOT be 1 or 2.
            if not has_edge and num_common_neighbors in [1, 2]:
                return False
    
    # If all pairs satisfy the conditions, the graph is a fixed point.
    return True

def find_fixed_point_graphs():
    """
    Finds and counts all non-isomorphic, connected graphs with at most 7
    vertices that are fixed points of the transformation T.
    """
    # Store counts for each size n
    counts_per_size = {}
    max_vertices = 7

    # networkx.graph_atlas_g() returns an iterator over all graphs in the
    # Atlas of Graphs, which includes all graphs up to 7 vertices.
    # This ensures we check one representative from each isomorphism class.
    all_graphs_in_atlas = nx.graph_atlas_g()

    for G in all_graphs_in_atlas:
        n = G.number_of_vertices()
        
        if n > max_vertices:
            continue
        
        # We only consider connected graphs as per the problem statement.
        if nx.is_connected(G):
            # Check if the graph is a fixed point of the transformation T.
            if check_graph_fixed_point(G):
                counts_per_size[n] = counts_per_size.get(n, 0) + 1

    # Format the final output as an equation, as requested.
    if not counts_per_size:
        print(0)
    else:
        # Sort the results by the number of vertices.
        sorted_counts = sorted(counts_per_size.items())
        
        equation_parts = [str(count) for n, count in sorted_counts]
        total_count = sum(counts_per_size.values())

        equation_str = " + ".join(equation_parts)
        
        # Only print the full equation if there's more than one term.
        if len(equation_parts) > 1:
            print(f"{equation_str} = {total_count}")
        else:
            print(total_count)

if __name__ == '__main__':
    find_fixed_point_graphs()