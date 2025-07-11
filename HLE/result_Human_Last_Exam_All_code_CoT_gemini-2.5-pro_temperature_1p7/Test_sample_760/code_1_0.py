import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs with at most 7 vertices
    that satisfy T(G) = G.
    """
    # graph_atlas_g() provides all non-isomorphic graphs up to 7 vertices
    graphs = nx.graph_atlas_g()
    
    found_graphs_count = 0
    # Use a dictionary to store counts for each vertex number n
    counts_per_n = {n: 0 for n in range(1, 8)}

    for G in graphs:
        # We only need to check graphs with at most 7 vertices, which the atlas provides.
        n = G.number_of_vertices()
        if n == 0:
            continue

        # Condition: The graph must be connected.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A of the graph G
        A = nx.to_numpy_array(G, dtype=int)
        
        # A single vertex graph with no edges.
        # A=[[0]], A_squared=[[0]]. No distinct pairs of vertices, so T(G)=G.
        if n == 1:
            found_graphs_count += 1
            counts_per_n[n] += 1
            continue
            
        # Calculate A_squared, which gives the number of length-2 paths
        A_squared = np.dot(A, A)
        
        # Create the adjacency matrix for the transformed graph T(G)
        A_transformed = np.zeros((n, n), dtype=int)
        
        for i in range(n):
            for j in range(i + 1, n):
                if A_squared[i, j] == 1 or A_squared[i, j] == 2:
                    A_transformed[i, j] = 1
                    A_transformed[j, i] = 1

        # Check if the transformation results in the same graph, i.e., T(G) = G
        if np.array_equal(A, A_transformed):
            found_graphs_count += 1
            counts_per_n[n] += 1

    print(f"The number of graphs found for each vertex count:")
    for n, count in counts_per_n.items():
        if count > 0:
            print(f"  - Vertices: {n}, Found: {count}")
            
    print("\nThe final calculation is:")
    
    # Collect the non-zero counts to form the equation string
    equation_numbers = [str(c) for c in counts_per_n.values() if c > 0]
    equation_str = " + ".join(equation_numbers)
    
    print(f"{equation_str} = {found_graphs_count}")

solve()