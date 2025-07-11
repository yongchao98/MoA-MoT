import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs with at most 7 vertices
    that satisfy T(G) = G, where T is the specified transformation.
    """

    # Plan:
    # 1. Iterate through all non-isomorphic graphs with up to 7 vertices.
    # 2. Filter for connected graphs.
    # 3. For each graph G, get its adjacency matrix A.
    # 4. Compute A^2 to find the number of length-2 paths.
    # 5. Construct the adjacency matrix A_T for the transformed graph T(G).
    # 6. Check if A_T is identical to A.
    # 7. Count the graphs that satisfy the condition.
    # 8. Print the counts for each size and the total.

    counts_by_size = {n: 0 for n in range(1, 8)}

    # nx.graph_atlas_g() provides all non-isomorphic graphs up to 7 vertices.
    for G in nx.graph_atlas_g():
        num_vertices = G.number_of_vertices()

        # Skip the null graph (0 vertices) and graphs with more than 7 vertices.
        if num_vertices == 0 or num_vertices > 7:
            continue

        # The problem specifies connected graphs.
        if not nx.is_connected(G):
            continue

        # For a single vertex graph, the condition is trivially met.
        if num_vertices == 1:
            counts_by_size[1] += 1
            continue
            
        # Get the adjacency matrix A of the graph G.
        # Sorting nodes ensures a consistent ordering.
        nodelist = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodelist)

        # Compute A^2. The entry (i, j) in A^2 is the number of length-2 paths
        # between vertex i and vertex j.
        A_squared = np.dot(A, A)

        # Construct the adjacency matrix A_T for the transformed graph T(G).
        A_T = np.zeros_like(A)
        for i in range(num_vertices):
            # Iterate over the upper triangle only, as the matrix is symmetric.
            for j in range(i + 1, num_vertices):
                num_paths = A_squared[i, j]
                # An edge exists in T(G) iff there are 1 or 2 length-2 paths.
                if num_paths == 1 or num_paths == 2:
                    A_T[i, j] = 1
                    A_T[j, i] = 1

        # Check if T(G) = G by comparing their adjacency matrices.
        if np.array_equal(A, A_T):
            counts_by_size[num_vertices] += 1

    # Output the results, showing each number in the final sum.
    equation_parts = []
    total_count = 0
    for n in range(1, 8):
        count = counts_by_size[n]
        print(f"Number of graphs with {n} vertices: {count}")
        total_count += count
        equation_parts.append(str(count))

    equation_str = " + ".join(equation_parts)
    print(f"Total number of graphs = {equation_str} = {total_count}")

solve()