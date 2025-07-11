import networkx as nx
import numpy as np

def solve_graph_problem():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 
    vertices that satisfy the condition T(G) = G.
    """
    # Initialize a dictionary to store the counts for each vertex size n.
    counts_by_size = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() returns a list of all non-isomorphic graphs
    # with up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    # Iterate through each graph to check the condition.
    for G in all_graphs:
        # The problem specifies connected graphs.
        if not nx.is_connected(G):
            continue
        
        n = G.number_of_nodes()
        if n == 0:
            continue
            
        # Get the adjacency matrix A.
        A = nx.to_numpy_array(G, dtype=int)
        
        # The number of length-2 paths is given by the matrix A^2.
        A_squared = np.dot(A, A)
        
        is_fixed_point = True
        # The condition T(G) = G means that for any two distinct vertices i and j:
        # an edge exists in G iff the number of length-2 paths between them is 1 or 2.
        # This translates to: A[i,j] == 1 <=> A_squared[i,j] in {1, 2}
        for i in range(n):
            for j in range(i + 1, n):
                edge_exists_in_G = (A[i, j] == 1)
                
                num_paths = A_squared[i, j]
                edge_exists_in_TG = (num_paths == 1 or num_paths == 2)
                
                if edge_exists_in_G != edge_exists_in_TG:
                    is_fixed_point = False
                    break
            if not is_fixed_point:
                break
                
        if is_fixed_point:
            counts_by_size[n] += 1
            
    # Output the results in the requested equation format.
    total_count = sum(counts_by_size.values())
    
    print("Number of valid graphs found for each vertex count n:")
    for n, count in counts_by_size.items():
        print(f"n = {n}: {count} graph(s)")
        
    equation_parts = [str(c) for c in counts_by_size.values()]
    equation_str = " + ".join(equation_parts)
    
    print("\nThe total number of graphs satisfying the condition is the sum:")
    print(f"{equation_str} = {total_count}")

solve_graph_problem()