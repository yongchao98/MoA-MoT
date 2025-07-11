import networkx as nx
import numpy as np

def solve_graph_problem():
    """
    Finds and counts all non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the condition T(G) = G.
    The transformation T connects two distinct vertices x and y if and only if
    there are exactly one or two paths of length 2 between them in G.
    """
    # Atlas of all non-isomorphic graphs with up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    # Dictionary to store the count of qualifying graphs for each vertex number.
    counts_per_n = {n: 0 for n in range(1, 8)}

    for G in all_graphs:
        n = len(G.nodes())

        # Skip the empty graph (n=0).
        if n == 0:
            continue

        # The problem requires the graph to be connected.
        if not nx.is_connected(G):
            continue

        # --- Check if T(G) = G ---

        # Case for n=1:
        # The graph has one vertex and no edges. It is connected.
        # The transformation rule applies to "distinct vertices x and y".
        # Since there are no such pairs, T(G) has no edges.
        # Thus, T(G) is the same as G.
        if n == 1:
            is_solution = True
        else:
            # For n >= 2, we use adjacency matrices.
            # Get the adjacency matrix A.
            A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))

            # Compute A^2, which counts paths of length 2.
            A_squared = A @ A

            # Construct the adjacency matrix A_T for the transformed graph T(G).
            A_T = np.zeros_like(A, dtype=int)
            for i in range(n):
                for j in range(i + 1, n):
                    num_paths_len_2 = A_squared[i, j]
                    if num_paths_len_2 == 1 or num_paths_len_2 == 2:
                        A_T[i, j] = 1
                        A_T[j, i] = 1
            
            # Check if T(G) = G by comparing their adjacency matrices.
            is_solution = np.array_equal(A, A_T)

        if is_solution:
            counts_per_n[n] += 1
            
    # --- Output the results ---
    print("Number of qualifying graphs found for each vertex count:")
    equation_parts = []
    for n in sorted(counts_per_n.keys()):
        count_at_n = counts_per_n[n]
        print(f"  - For {n} vertices: {count_at_n} graph(s)")
        equation_parts.append(str(count_at_n))
        
    total_count = sum(counts_per_n.values())
    
    print("\nThe final count is the sum of these numbers:")
    final_equation = " + ".join(equation_parts)
    print(f"{final_equation} = {total_count}")

# Execute the solution
solve_graph_problem()
<<<4>>>