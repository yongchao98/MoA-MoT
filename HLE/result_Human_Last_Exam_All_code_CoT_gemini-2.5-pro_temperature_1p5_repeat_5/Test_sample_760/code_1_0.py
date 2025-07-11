import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the transformation property T(G) = G.
    """
    
    # We will store the counts for each number of vertices n.
    counts_by_n = {n: 0 for n in range(1, 8)}

    # The networkx.graph_atlas_g() function provides a list of all non-isomorphic
    # graphs on 0 to 7 vertices.
    graph_list = nx.graph_atlas_g()

    for G in graph_list:
        n = G.number_of_vertices()

        # Skip the empty graph (n=0) or graphs larger than 7 vertices.
        if n == 0 or n > 7:
            continue

        # The problem specifies that the graphs must be connected.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A of the graph G.
        A = nx.to_numpy_array(G, dtype=int)
        
        # Compute A^2. The entry (A^2)_ij is the number of length-2 paths
        # between vertices i and j.
        A_sq = np.dot(A, A)

        # Build the adjacency matrix A_T for the transformed graph T(G).
        A_T = np.zeros_like(A)
        for i in range(n):
            for j in range(i + 1, n):
                # The rule: add an edge if there are 1 or 2 length-2 paths.
                if A_sq[i, j] == 1 or A_sq[i, j] == 2:
                    A_T[i, j] = 1
                    A_T[j, i] = 1

        # Check if the transformation results in the same graph, i.e., T(G) = G.
        # This is true if their adjacency matrices are identical.
        if np.array_equal(A, A_T):
            counts_by_n[n] += 1
    
    total_count = 0
    equation_parts = []
    
    for n in range(1, 8):
        count = counts_by_n[n]
        print(f"Found {count} graph(s) with {n} vertices.")
        total_count += count
        equation_parts.append(str(count))

    equation_str = " + ".join(equation_parts)
    print("\nFinal calculation:")
    print(f"{equation_str} = {total_count}")


solve()