import networkx as nx
import numpy as np

def find_fixed_point_graphs():
    """
    Finds the number of non-isomorphic, connected graphs with at most 7 vertices
    that are fixed points of the transformation T.

    The transformation T is defined as follows: for two distinct vertices x and y,
    an edge is added between x and y in T(G) if and only if there are exactly
    one or two length-2 paths between x and y in G.
    """
    total_count = 0
    
    print("Counting graphs G such that T(G) = G for n vertices...")
    print("-------------------------------------------------------")

    # Iterate through the number of vertices from 1 to 7
    for n in range(1, 8):
        count_for_n = 0
        
        # Generate all non-isomorphic graphs with n vertices
        graph_generator = nx.nonisomorphic_graphs(n)
        
        for G in graph_generator:
            # We are only interested in connected graphs
            if not nx.is_connected(G):
                continue
            
            # Get the adjacency matrix A of the graph G
            # Using a sorted nodelist ensures a canonical representation
            nodelist = sorted(G.nodes())
            A = nx.to_numpy_array(G, nodelist=nodelist)
            
            # The number of length-2 paths between vertices is given by A^2
            A_squared = np.dot(A, A)
            
            # Construct the adjacency matrix A_prime for the transformed graph T(G)
            # An edge (i, j) exists in T(G) if the number of length-2 paths
            # is 1 or 2.
            A_prime = ((A_squared >= 1) & (A_squared <= 2)).astype(int)
            
            # The transformation does not create self-loops, so the diagonal must be zero.
            np.fill_diagonal(A_prime, 0)
            
            # Check if the transformation is a fixed point (i.e., T(G) = G)
            # This is equivalent to their adjacency matrices being equal.
            if np.array_equal(A, A_prime):
                count_for_n += 1

        print(f"Number of graphs with {n} vertices: {count_for_n}")
        total_count += count_for_n
        
    print("-------------------------------------------------------")
    print(f"The total number of such graphs is the sum: {total_count}")

if __name__ == '__main__':
    find_fixed_point_graphs()