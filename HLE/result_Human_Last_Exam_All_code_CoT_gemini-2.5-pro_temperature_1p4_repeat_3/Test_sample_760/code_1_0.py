import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 
    vertices that are fixed points of the transformation T, i.e., T(G) = G.
    """
    # Store the counts for each number of vertices n
    counts_per_size = {n: 0 for n in range(1, 8)}
    
    # Iterate through the number of vertices from 1 to 7
    for n in range(1, 8):
        # Generate all non-isomorphic graphs of size n
        graph_generator = nx.nonisomorphic_graphs(n)
        
        for G in graph_generator:
            # The problem states the graph must be connected
            if not nx.is_connected(G):
                continue
            
            # The graph with one vertex (n=1) has no edges.
            # A is [[0]], A^2 is [[0]]. The condition is vacuously true.
            if n == 1:
                counts_per_size[n] += 1
                continue

            # Get the adjacency matrix A for the graph G
            nodelist = sorted(G.nodes())
            A = nx.to_numpy_array(G, nodelist=nodelist)
            
            # Compute A^2, which gives the number of 2-paths
            A_squared = np.dot(A, A)
            
            is_fixed_point = True
            # Check the condition T(G) = G for all pairs of distinct vertices
            for i in range(n):
                for j in range(i + 1, n):
                    # Check if an edge exists between i and j in the original graph G
                    edge_exists_in_G = (A[i, j] == 1)
                    
                    # Count the number of 2-paths between i and j in G
                    num_2_paths = A_squared[i, j]
                    
                    # Determine if an edge would exist in the transformed graph T(G)
                    edge_exists_in_TG = (1 <= num_2_paths <= 2)
                    
                    # For T(G) = G, the existence of an edge must be the same
                    if edge_exists_in_G != edge_exists_in_TG:
                        is_fixed_point = False
                        break
                if not is_fixed_point:
                    break
            
            if is_fixed_point:
                counts_per_size[n] += 1
    
    # Format and print the final equation and result
    equation_parts = [str(counts_per_size[n]) for n in range(1, 8)]
    total_count = sum(counts_per_size.values())
    equation_str = " + ".join(equation_parts) + f" = {total_count}"
    
    print(equation_str)

solve()