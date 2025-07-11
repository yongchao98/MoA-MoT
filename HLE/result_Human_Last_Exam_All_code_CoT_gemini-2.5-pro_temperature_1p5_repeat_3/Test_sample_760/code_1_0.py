import networkx as nx
import numpy as np

def solve_graph_transformation():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the condition T(G) = G.

    The transformation T(G) creates an edge between vertices x and y iff there are
    1 or 2 paths of length 2 between them in G.

    The condition T(G) = G is checked by comparing the graph's adjacency matrix A
    with the matrix derived from A^2. Specifically, for any distinct vertices i and j,
    A[i,j] must be 1 if and only if A^2[i,j] is 1 or 2.
    """

    def is_fixed_point(G):
        """
        Checks if a graph G is a fixed point of the transformation T.
        """
        n = G.number_of_vertices()

        # A graph with 0 or 1 vertex has no pairs of distinct vertices, so the
        # condition is vacuously true. The transformed graph T(G) has an empty
        # edge set, same as G. So it's a solution.
        if n <= 1:
            return True

        # Use a sorted nodelist to ensure consistent matrix representation.
        nodelist = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodelist, dtype=int)
        
        # Compute the square of the adjacency matrix. A_squared[i, j] counts
        # the number of paths of length 2 between i and j.
        A_squared = np.dot(A, A)

        # Check the condition for all pairs of distinct vertices (i, j).
        for i in range(n):
            for j in range(i + 1, n):
                has_edge_in_G = (A[i, j] == 1)
                
                num_paths = A_squared[i, j]
                has_edge_in_TG = (num_paths == 1 or num_paths == 2)
                
                # The condition T(G) = G fails if the edge status differs.
                if has_edge_in_G != has_edge_in_TG:
                    return False
        
        return True

    # Store the count of solutions for each number of vertices.
    counts_per_n = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() provides all non-isomorphic graphs up to n=7.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_vertices()

        if n > 7:
            break
        if n == 0:
            continue
        
        # The problem requires connected graphs.
        if nx.is_connected(G):
            if is_fixed_point(G):
                counts_per_n[n] += 1

    total_count = sum(counts_per_n.values())

    # Build and print the equation as requested.
    equation_parts = [f"{counts_per_n[n]}" for n in range(1, 8)]
    equation_str = " + ".join(equation_parts)

    print(f"The number of such graphs for n=1 to 7 are, respectively:")
    print(f"{equation_str}")
    print(f"Total number of graphs = {total_count}")

solve_graph_transformation()
<<<3>>>