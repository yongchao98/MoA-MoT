import networkx as nx
import numpy as np

def solve_graph_transformation_problem():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that are fixed points of the transformation T, i.e., T(G) = G.
    """
    # A dictionary to store the count of solutions for each number of vertices.
    counts_by_n = {i: 0 for i in range(1, 8)}

    # networkx.graph_atlas_g() provides a list of all non-isomorphic graphs
    # with up to 7 vertices.
    graph_atlas = nx.graph_atlas_g()

    for G in graph_atlas:
        n = G.number_of_nodes()

        # Skip the null graph (0 vertices) and graphs with more than 7 vertices.
        if n == 0 or n > 7:
            continue

        # The problem specifies connected graphs.
        if not nx.is_connected(G):
            continue

        # Handle the n=1 case separately. A graph with one vertex is connected.
        # The transformation T is vacuously satisfied as there are no pairs of
        # distinct vertices. T(G) has no edges, just like G. So, it's a solution.
        if n == 1:
            counts_by_n[n] += 1
            continue

        # Get the adjacency matrix A of the graph G.
        A = nx.to_numpy_array(G, dtype=int)

        # Compute A^2. The entry (A^2)_ij is the number of length-2 paths from i to j.
        A_squared = np.dot(A, A)

        # Construct the adjacency matrix A_prime for the transformed graph T(G).
        # We can use boolean masking for a concise implementation.
        # First, ensure diagonal entries of A_squared won't fall in the [1, 2] range.
        np.fill_diagonal(A_squared, 0)
        
        # An edge exists in T(G) if the number of length-2 paths is 1 or 2.
        A_prime = ((A_squared >= 1) & (A_squared <= 2)).astype(int)

        # Check if G is a fixed point by comparing its adjacency matrix with T(G)'s.
        if np.array_equal(A, A_prime):
            counts_by_n[n] += 1

    # Calculate the total count.
    total_count = sum(counts_by_n.values())

    # Print the results, including the breakdown as requested.
    print("The number of solutions for each number of vertices is:")
    for n, count in counts_by_n.items():
        print(f"n = {n}: {count} graph(s)")

    print("\nThe final equation is:")
    equation_str = " + ".join(str(c) for c in counts_by_n.values())
    print(f"{equation_str} = {total_count}")

solve_graph_transformation_problem()