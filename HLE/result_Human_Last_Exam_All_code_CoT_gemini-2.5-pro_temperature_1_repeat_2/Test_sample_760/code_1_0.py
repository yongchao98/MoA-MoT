import networkx as nx
import numpy as np

def check_graph_is_solution(G):
    """
    Checks if a graph G satisfies the condition T(G) = G.
    """
    n = G.number_of_vertices()
    
    # A graph with 0 vertices is not considered.
    # A graph with 1 vertex has no pairs of distinct vertices, so the condition
    # is vacuously true. It is connected.
    if n == 0:
        return False
    if n == 1:
        return True
    
    # Get the adjacency matrix A
    A = nx.to_numpy_array(G, dtype=int)
    
    # Compute A^2
    A_squared = np.dot(A, A)
    
    # Check the condition for all distinct pairs of vertices (i, j)
    for i in range(n):
        for j in range(i + 1, n):
            # Check if an edge exists in the original graph G
            has_edge = (A[i, j] == 1)
            
            # Get the number of length-2 paths
            num_paths = A_squared[i, j]
            
            # Check if an edge would exist in the transformed graph T(G)
            transformed_has_edge = (num_paths == 1 or num_paths == 2)
            
            # For T(G) = G, the edge existence must be identical
            if has_edge != transformed_has_edge:
                return False
                
    return True

def solve():
    """
    Finds and counts the number of non-isomorphic, connected graphs G 
    with at most 7 vertices such that T(G) = G.
    """
    # Store the counts for each number of vertices n
    counts_by_n = {n: 0 for n in range(1, 8)}
    
    # networkx.graph_atlas_g() returns a list of all non-isomorphic graphs
    # on up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_vertices()
        
        # We only consider graphs with at most 7 vertices
        if n > 7:
            continue
            
        # The problem specifies connected graphs
        if not nx.is_connected(G):
            continue

        # Check if the graph is a solution
        if check_graph_is_solution(G):
            counts_by_n[n] += 1
            
    # Output the results as requested
    total_count = 0
    equation_parts = []
    
    print("Finding the number of graphs G such that T(G) = G:")
    for n in sorted(counts_by_n.keys()):
        count = counts_by_n[n]
        if count > 0:
            print(f"Number of graphs with {n} vertices: {count}")
            total_count += count
            equation_parts.append(str(count))
    
    equation_str = " + ".join(equation_parts)
    print(f"\nTotal number of graphs = {equation_str} = {total_count}")

# Execute the solver
solve()