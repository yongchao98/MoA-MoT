import networkx as nx
import numpy as np

def solve_graph_problem():
    """
    This function solves the graph theory problem by iterating through all
    non-isomorphic, connected graphs with 1 to 7 vertices. It checks each
    graph to see if it is a fixed point of the transformation T.
    """
    counts_per_n = []
    
    # Iterate through the number of vertices from 1 to 7
    for n in range(1, 8):
        n_count = 0
        # Generate all non-isomorphic, connected graphs of size n
        graphs_iterator = nx.connected_graphs_iterator(n)
        
        for g in graphs_iterator:
            is_solution = False
            # The single-vertex graph is a trivial solution.
            if n == 1:
                is_solution = True
            else:
                # Get the adjacency matrix A
                nodes = sorted(g.nodes())
                A = nx.to_numpy_array(g, nodelist=nodes)
                
                # Compute A^2, which gives the number of length-2 paths
                A2 = np.dot(A, A)
                
                # Construct the adjacency matrix A_prime for the transformed graph T(G)
                A_prime = np.zeros_like(A)
                
                # An edge exists in T(G) if the number of length-2 paths is 1 or 2.
                condition_met = (A2 >= 1) & (A2 <= 2)
                A_prime[condition_met] = 1
                
                # Ensure no self-loops in the transformed graph
                np.fill_diagonal(A_prime, 0)
                
                # Check if the original graph G is a fixed point of T (i.e., A = A')
                if np.array_equal(A, A_prime):
                    is_solution = True

            if is_solution:
                n_count += 1
        
        counts_per_n.append(n_count)

    # Print the results in the required format
    total_count = sum(counts_per_n)
    
    print("Number of satisfying graphs found for each number of vertices:")
    for i, count in enumerate(counts_per_n):
        num_vertices = i + 1
        print(f"n = {num_vertices}: {count}")
        
    print("\nThe final calculation is:")
    equation_str = " + ".join(map(str, counts_per_n))
    print(f"{equation_str} = {total_count}")

solve_graph_problem()