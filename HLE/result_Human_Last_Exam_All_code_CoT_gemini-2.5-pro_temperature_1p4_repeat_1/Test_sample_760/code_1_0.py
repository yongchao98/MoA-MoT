import networkx as nx
import numpy as np

def solve_graph_problem():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy T(G) = G, where T is the transformation defined in the problem.
    """
    # Dictionary to store the count of fixed-point graphs for each number of vertices
    counts_per_n = {n: 0 for n in range(1, 8)}
    
    # nx.graph_atlas_g() returns a list of all non-isomorphic graphs on up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_vertices()
        
        # We are interested in graphs with 1 to 7 vertices.
        if n == 0 or n > 7:
            continue
        
        # The problem specifies that the graphs must be connected.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A. Ensure consistent node ordering for matrix operations.
        node_list = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=node_list)
        
        # Compute A^2. The entry (i, j) of A^2 gives the number of length-2 paths between i and j.
        A2 = np.dot(A, A)
        
        is_fixed_point = True
        # For n=1, the loops won't run, and is_fixed_point remains True, which is correct.
        for i in range(n):
            for j in range(i + 1, n):
                # Check if there is an edge in the original graph G.
                has_edge = (A[i, j] == 1)
                
                # Get the number of length-2 paths between i and j.
                num_paths_len_2 = A2[i, j]
                
                # According to the transformation T, an edge should exist in T(G) if there are 1 or 2 paths of length 2.
                t_has_edge = (1 <= num_paths_len_2 <= 2)
                
                # For G to be a fixed point (T(G) = G), the edge existence in G must match the condition for T(G).
                if has_edge != t_has_edge:
                    is_fixed_point = False
                    break
            if not is_fixed_point:
                break
        
        if is_fixed_point:
            counts_per_n[n] += 1
            
    total_count = sum(counts_per_n.values())
    
    print("Number of graphs found for each vertex count 'n':")
    for n in sorted(counts_per_n.keys()):
        print(f"n = {n}: {counts_per_n[n]} graphs")
    
    print("\nThe total count is the sum of these numbers.")
    # Output each number in the final equation, as requested.
    sum_str_parts = [str(counts_per_n[n]) for n in sorted(counts_per_n.keys())]
    sum_str = " + ".join(sum_str_parts)
    print(f"Final Calculation: {sum_str} = {total_count}")

if __name__ == "__main__":
    solve_graph_problem()