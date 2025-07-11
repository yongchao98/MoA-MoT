import networkx as nx
import numpy as np

def solve_graph_problem():
    """
    Finds and counts the number of non-isomorphic, connected graphs G with at most 7 vertices
    that are fixed points of the transformation T, i.e., T(G) = G.
    """
    # Array to store the count of solutions for each number of vertices n=1..7
    solution_counts_by_n = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() provides all non-isomorphic graphs up to 7 vertices.
    # The graphs are sorted by number of nodes.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        # We only consider graphs with 1 to 7 vertices.
        if n == 0:
            continue
        if n > 7:
            break

        # The problem requires the graphs to be connected.
        if not nx.is_connected(G):
            continue

        # Case n=1: A single vertex graph. It's connected.
        # The transformation T is defined for "distinct vertices x and y".
        # Since there are no such pairs, the condition for adding edges is never met.
        # T(G) has no edges, just like G. So T(G) = G.
        if n == 1:
            solution_counts_by_n[1] += 1
            continue

        # For n > 1, we perform the matrix-based check.
        # Get the adjacency matrix A. A sorted nodelist ensures consistent ordering.
        nodelist = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodelist)

        # Compute A^2. The entry (i, j) is the number of length-2 paths.
        A_squared = np.dot(A, A)

        # Build the adjacency matrix A_T for the transformed graph T(G).
        A_T = np.zeros_like(A, dtype=int)
        
        for i in range(n):
            for j in range(i + 1, n):
                # An edge exists in T(G) iff there are 1 or 2 length-2 paths.
                if A_squared[i, j] == 1 or A_squared[i, j] == 2:
                    A_T[i, j] = 1
                    A_T[j, i] = 1
        
        # Check if T(G) = G by comparing their adjacency matrices.
        if np.array_equal(A, A_T):
            solution_counts_by_n[n] += 1

    # Print the results in a clear format.
    print("Number of graphs G satisfying T(G) = G, for each number of vertices n:")
    total_solutions = 0
    equation_parts = []
    for n in range(1, 8):
        count = solution_counts_by_n[n]
        print(f"n = {n}: {count} graphs")
        if count > 0:
            total_solutions += count
            # This part is for the final equation format.
            # Showing the actual numbers that sum up.
            equation_parts.append(str(count))

    print("\nThe total number of such graphs is the sum of the counts for each n:")
    # We output each number in the final equation as requested.
    equation_str = " + ".join(f"{solution_counts_by_n[n]}" for n in range(1, 8))
    print(f"{equation_str} = {total_solutions}")

# Run the solver.
solve_graph_problem()