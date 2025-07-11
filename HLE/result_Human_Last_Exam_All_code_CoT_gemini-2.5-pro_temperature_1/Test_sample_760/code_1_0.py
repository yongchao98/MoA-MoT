import networkx as nx
import numpy as np

def find_stable_graphs():
    """
    This function finds the number of non-isomorphic, connected graphs G with at most 7
    vertices that are stable under the transformation T, i.e., T(G) = G.

    The transformation T is defined as follows: for two distinct vertices x and y,
    an edge is added between x and y in T(G) if and only if there are exactly
    one or two length-2 paths between x and y in G.
    """
    
    # A dictionary to store the count of solution graphs for each number of vertices n.
    solution_counts_by_n = {i: 0 for i in range(1, 8)}

    # networkx.graph_atlas_g() provides a list of all non-isomorphic graphs with up to 7 vertices.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        # Skip the empty graph (n=0)
        if n == 0:
            continue

        # The problem requires the graphs to be connected.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A of the graph G.
        # Ensure the data type is suitable for matrix multiplication.
        A = nx.to_numpy_array(G, dtype=int)

        # Calculate A^2. The entry (A^2)[i, j] is the number of length-2 paths
        # between vertex i and vertex j.
        A_sq = np.dot(A, A)
        
        # Construct the adjacency matrix A_T for the transformed graph T(G).
        # An edge (i, j) exists in T(G) if 1 <= (A^2)[i, j] <= 2.
        # We can use boolean masking and type casting to create A_T efficiently.
        A_T = ((A_sq >= 1) & (A_sq <= 2)).astype(int)
        
        # The transformation is defined for distinct vertices, which means T(G)
        # has no self-loops. The diagonal of its adjacency matrix must be all zeros.
        # A_sq[i, i] is the degree of vertex i, which can be 1 or 2. The boolean
        # mask would create 1s on the diagonal, so we must reset them to 0.
        np.fill_diagonal(A_T, 0)
        
        # Check if the transformation results in the same graph, i.e., T(G) = G.
        # This is true if their adjacency matrices are identical.
        if np.array_equal(A, A_T):
            solution_counts_by_n[n] += 1

    # Print the results in a clear format.
    print("Finding the number of graphs G where T(G) = G for n <= 7.")
    print("-" * 50)
    for n, count in solution_counts_by_n.items():
        print(f"Number of solutions for n = {n} vertices: {count}")
    print("-" * 50)

    # The final output is the sum of counts for n=1 to 7.
    # The instruction "output each number in the final equation" suggests this format.
    counts = list(solution_counts_by_n.values())
    total_solutions = sum(counts)
    
    equation_str = " + ".join(map(str, counts))
    
    print("The total count is the sum of counts for each number of vertices:")
    print(f"{equation_str} = {total_solutions}")

if __name__ == '__main__':
    find_stable_graphs()