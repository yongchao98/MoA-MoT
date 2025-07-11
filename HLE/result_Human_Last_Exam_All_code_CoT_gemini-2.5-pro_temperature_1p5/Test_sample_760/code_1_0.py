import networkx as nx
import numpy as np

def find_fixed_point_graphs():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7
    vertices that are fixed points of a specific transformation T.

    The transformation T(G) has an edge between vertices x and y iff there are
    1 or 2 paths of length 2 between x and y in G.
    """
    counts_per_n = []
    
    # Iterate through the number of vertices from 1 to 7
    for n in range(1, 8):
        count_for_n = 0
        # networkx.graph_atlas(n) returns a list of all non-isomorphic graphs with n vertices.
        graphs = nx.graph_atlas(n)
        
        for g in graphs:
            # Skip graphs that are not connected
            if not nx.is_connected(g):
                continue

            # Get the adjacency matrix A of the graph G
            # Using a sorted nodelist ensures consistent matrix representation
            A = nx.to_numpy_array(g, nodelist=sorted(g.nodes()))

            # The number of length-2 paths between vertices is given by A^2
            A_squared = np.dot(A, A)

            # Construct the adjacency matrix A_T for the transformed graph T(G)
            # An edge (i,j) exists in T(G) if A_squared[i,j] is 1 or 2.
            A_T = ((A_squared >= 1) & (A_squared <= 2)).astype(int)
            
            # The transformation does not create self-loops, so set the diagonal to 0
            np.fill_diagonal(A_T, 0)

            # Check if the transformation results in the original graph (T(G) = G)
            if np.array_equal(A, A_T):
                count_for_n += 1
        
        counts_per_n.append(count_for_n)

    # Print the results in the format: N1 + N2 + ... + N7 = Total
    equation_str = " + ".join(map(str, counts_per_n))
    total_count = sum(counts_per_n)
    print(f"{equation_str} = {total_count}")

# Execute the function to find and print the result
find_fixed_point_graphs()