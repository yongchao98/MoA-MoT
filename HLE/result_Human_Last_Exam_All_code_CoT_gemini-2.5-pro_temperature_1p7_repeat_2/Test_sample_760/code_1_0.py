import networkx as nx
import numpy as np

def solve_graph_transformation():
    """
    Finds and counts the number of non-isomorphic, connected graphs G
    with at most 7 vertices that satisfy the condition T(G) = G.
    """
    # This dictionary will store the count of qualifying graphs for each number of vertices.
    counts_by_n = {i: 0 for i in range(1, 8)}

    # The networkx graph atlas contains all non-isomorphic graphs up to 7 vertices.
    # We iterate through them to find the ones that meet the criteria.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_nodes()

        # Skip the empty graph (n=0). The atlas is capped at 7 vertices.
        if n == 0:
            continue

        # The problem specifies using connected graphs.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix 'A' of the graph G.
        A = nx.to_numpy_array(G, dtype=np.int32)

        # The number of paths of length 2 between any two vertices is given by A^2.
        A2 = np.dot(A, A)

        # Construct the adjacency matrix 'A_T' for the transformed graph T(G).
        # An edge exists in T(G) if there are 1 or 2 length-2 paths in G.
        A_T = np.zeros_like(A)
        # Iterate over the upper triangle of the matrix as it's symmetric.
        for i in range(n):
            for j in range(i + 1, n):
                if A2[i, j] == 1 or A2[i, j] == 2:
                    A_T[i, j] = 1
                    A_T[j, i] = 1

        # Check if T(G) = G by comparing their adjacency matrices.
        if np.array_equal(A, A_T):
            counts_by_n[n] += 1

    # Print the results in a structured way, showing the breakdown by number of vertices.
    total_count = 0
    equation_parts = []
    print("Finding the number of graphs G for each vertex count n where T(G) = G:")
    for n in sorted(counts_by_n.keys()):
        count = counts_by_n[n]
        print(f"Number of graphs with {n} vertices: {count}")
        total_count += count
        equation_parts.append(str(count))

    # To satisfy the "output each number in the final equation" hint, we show the sum.
    equation_str = " + ".join(equation_parts)
    print("\nThe final calculation is:")
    print(f"{equation_str} = {total_count}")


solve_graph_transformation()
<<<5>>>