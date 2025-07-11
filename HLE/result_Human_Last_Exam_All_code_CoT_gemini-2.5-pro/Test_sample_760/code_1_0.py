import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 
    vertices that satisfy the condition T(G) = G.
    """
    
    total_count = 0
    # A dictionary to store the count of valid graphs for each number of vertices.
    counts_per_n = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() provides a list of all non-isomorphic graphs
    # with up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_nodes()

        # We skip the null graph (0 vertices) and only consider connected graphs.
        if n > 0 and nx.is_connected(G):
            # Get the adjacency matrix A of the graph G.
            A = nx.to_numpy_array(G, dtype=int)

            # Compute A^2. The (i, j)-th entry gives the number of length-2 paths.
            A_squared = np.dot(A, A)
            
            is_match = True
            # Check if the condition T(G) = G holds for all pairs of vertices.
            for i in range(n):
                for j in range(i + 1, n):
                    num_common_neighbors = A_squared[i, j]
                    
                    # An edge exists in the transformed graph T(G) if there are 
                    # 1 or 2 common neighbors (length-2 paths).
                    edge_in_T_G = (num_common_neighbors == 1 or num_common_neighbors == 2)
                    
                    # An edge exists in the original graph G.
                    edge_in_G = (A[i, j] == 1)
                    
                    # If the edge status in G and T(G) is different, the graph
                    # does not satisfy the condition.
                    if edge_in_G != edge_in_T_G:
                        is_match = False
                        break
                if not is_match:
                    break
            
            if is_match:
                counts_per_n[n] += 1

    # Print the results in the requested format.
    print("The number of graphs G satisfying T(G)=G for n vertices:")
    equation_parts = []
    for n in sorted(counts_per_n.keys()):
        count_n = counts_per_n[n]
        print(f"For n = {n}, there are {count_n} graphs.")
        total_count += count_n
        equation_parts.append(str(count_n))

    print("\nThe final equation for the total count is:")
    final_equation = " + ".join(equation_parts) + f" = {total_count}"
    print(final_equation)

solve()