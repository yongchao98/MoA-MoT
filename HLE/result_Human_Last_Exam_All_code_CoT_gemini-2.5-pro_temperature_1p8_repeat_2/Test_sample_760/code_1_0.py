import networkx as nx
import numpy as np

def find_fixed_point_graphs():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that are fixed points of the transformation T, i.e., T(G) = G.

    The transformation T is defined as: an edge exists between distinct vertices
    x and y in T(G) iff there are exactly one or two length-2 paths between
    x and y in G.
    """
    total_count = 0
    # Store the counts for each number of vertices
    counts_per_n = {}

    # Iterate through the number of vertices from 1 to 7
    for n in range(1, 8):
        count_for_n = 0
        # Using a generator for non-isomorphic graphs is efficient.
        # It's guaranteed to work for n <= 8 in NetworkX.
        graph_generator = nx.nonisomorphic_graphs(n)

        for G in graph_generator:
            # The problem specifies connected graphs only
            if not nx.is_connected(G):
                continue

            # Get the adjacency matrix A.
            # Using a sorted nodelist ensures a canonical matrix representation.
            nodes = sorted(G.nodes())
            A = nx.to_numpy_array(G, nodelist=nodes)

            # Compute A^2. The entry (i, j) gives the number of length-2 paths.
            A_squared = np.dot(A, A)

            # Construct the adjacency matrix A_prime for the transformed graph T(G)
            # The condition is that the number of length-2 paths is 1 or 2.
            condition = (A_squared == 1) | (A_squared == 2)
            A_prime = condition.astype(int)
            
            # The transformation is defined for distinct vertices, so no self-loops.
            # We set the diagonal of A_prime to zero.
            np.fill_diagonal(A_prime, 0)
            
            # Check if the transformation results in the same graph (T(G) = G)
            if np.array_equal(A, A_prime):
                count_for_n += 1
        
        if count_for_n > 0:
            counts_per_n[n] = count_for_n
        total_count += count_for_n
    
    # --- Output the results as requested ---
    print("Found graphs satisfying the condition for the following vertex counts:")
    equation_parts = []
    # Sort by n for a clean output
    for n_val, count_val in sorted(counts_per_n.items()):
        print(f"Number of vertices (n={n_val}): {count_val} graph(s)")
        equation_parts.append(str(count_val))

    # To satisfy the "output each number in the final equation" instruction
    final_equation = " + ".join(equation_parts) + f" = {total_count}"
    print("\nThe final count is the sum of these numbers:")
    print(final_equation)

# Execute the function to find the answer
find_fixed_point_graphs()