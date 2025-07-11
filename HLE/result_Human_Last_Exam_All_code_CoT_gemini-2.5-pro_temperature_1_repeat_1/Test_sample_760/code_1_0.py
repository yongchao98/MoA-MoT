import networkx as nx
import numpy as np

def find_fixed_point_graphs():
    """
    Finds and counts the number of non-isomorphic, connected graphs G with at most 7
    vertices that are fixed points of the transformation T, i.e., T(G) = G.
    
    The transformation T is defined as follows: an edge exists between distinct 
    vertices x and y in T(G) if and only if there are exactly one or two 
    length-2 paths between x and y in G.
    """
    
    # Store the count of solutions for each number of vertices n.
    counts_per_n = {n: 0 for n in range(1, 8)}

    # Handle the n=1 case (graph with a single vertex).
    # It is connected. Its adjacency matrix is A = [[0]].
    # A^2 is also [[0]]. The condition for adding edges applies to distinct vertices,
    # so the transformed graph also has a [[0]] adjacency matrix. Thus, T(G) = G.
    counts_per_n[1] = 1

    # nx.graph_atlas_g() returns a list of all non-isomorphic graphs up to 7 vertices.
    # It starts with the null graph (0 vertices), so we filter for n in [2, 7].
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_vertices()

        if not (2 <= n <= 7):
            continue

        # We only consider connected graphs.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A of the graph G.
        # Ensure a consistent node ordering for matrix comparison.
        nodes = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodes)

        # The number of length-2 paths is given by A^2.
        A_sq = np.dot(A, A)

        # Construct the adjacency matrix A_T for the transformed graph T(G).
        # An edge exists if the number of length-2 paths is 1 or 2.
        A_T = ((A_sq >= 1) & (A_sq <= 2)).astype(int)

        # Set diagonal to zero, as self-loops are not allowed.
        np.fill_diagonal(A_T, 0)

        # Check if the transformed graph is identical to the original graph.
        if np.array_equal(A, A_T):
            counts_per_n[n] += 1

    # Output the results in the required format.
    print("Number of graphs G satisfying T(G) = G for each number of vertices n:")
    total_count = 0
    equation_parts = []
    for n in range(1, 8):
        count = counts_per_n[n]
        print(f"For n = {n}, found {count} graph(s).")
        total_count += count
        equation_parts.append(str(count))

    print("\nThe total number of such graphs is the sum of the counts for each n:")
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_count}")

    # Final answer in the specified format
    print(f"\n<<<{total_count}>>>")

if __name__ == '__main__':
    find_fixed_point_graphs()