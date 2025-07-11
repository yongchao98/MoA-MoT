import networkx as nx
import numpy as np

def find_fixed_point_graphs():
    """
    Finds and counts non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the condition T(G) = G.

    The transformation T is defined as follows: for two distinct vertices x and y,
    an edge is added between x and y in T(G) if and only if there are exactly
    one or two length-2 paths between x and y in G.
    """
    
    # Store the count of solutions for each number of vertices n.
    counts_by_n = {i: 0 for i in range(1, 8)}

    # nx.graph_atlas_g() returns a list of all non-isomorphic graphs with up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_vertices()

        # We are interested in graphs with at most 7 vertices. The atlas stops at 7.
        # We also need to skip the null graph (n=0).
        if n == 0:
            continue

        # The problem requires the graphs to be connected.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A of the graph G.
        A = nx.to_numpy_array(G, dtype=int)

        # The number of length-2 paths between vertices is given by A^2.
        A_squared = np.dot(A, A)

        # Construct the adjacency matrix for the transformed graph T(G).
        # Start with a zero matrix.
        T_A = np.zeros_like(A)

        # An edge (i, j) exists in T(G) if 1 <= (A^2)[i, j] <= 2.
        condition_mask = (A_squared >= 1) & (A_squared <= 2)
        T_A[condition_mask] = 1

        # The transformation is for distinct vertices, so no self-loops.
        # This step also correctly handles the diagonal of A^2, which counts degrees.
        np.fill_diagonal(T_A, 0)

        # Check if the transformation results in the original graph (T(G) = G).
        if np.array_equal(A, T_A):
            counts_by_n[n] += 1
            
    # Format and print the final equation showing the breakdown by vertex count.
    total_count = sum(counts_by_n.values())
    equation_parts = [str(counts_by_n[i]) for i in range(1, 8)]
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_count}")

# Execute the function to find the graphs and print the result.
find_fixed_point_graphs()