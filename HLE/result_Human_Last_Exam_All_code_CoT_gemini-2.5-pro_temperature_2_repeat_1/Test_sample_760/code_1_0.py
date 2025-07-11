import networkx as nx
import numpy as np

def solve_graph_transformation():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy T(G) = G, where T is a specific graph transformation.

    An edge (x, y) exists in T(G) if and only if there are exactly one or two
    length-2 paths between x and y in G.
    """
    
    # This dictionary will store the counts of found graphs, keyed by the number of vertices.
    counts_by_n = {}
    
    # networkx.graph_atlas_g() returns a list of all non-isomorphic graphs
    # with up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_vertices()

        # We are looking for connected graphs with at least one vertex.
        if n == 0 or not nx.is_connected(G):
            continue

        # Get the adjacency matrix of the graph G.
        A = nx.to_numpy_array(G, dtype=int)
        
        # The number of length-2 paths (common neighbors) between any two vertices i and j
        # is given by the (i, j)-th entry of the matrix A^2.
        A2 = np.matmul(A, A)
        
        # Construct the adjacency matrix A_prime for the transformed graph T(G).
        # An edge (i, j) exists in T(G) iff A2[i, j] is 1 or 2.
        is_one_or_two = (A2 == 1) | (A2 == 2)
        A_prime = is_one_or_two.astype(int)
        
        # The transformation T does not create self-loops, so the diagonal of A_prime must be zero.
        # Our source graph G has no self-loops (A has zero diagonal). A^2 has vertex degrees
        # on its diagonal, so the boolean mask could make diagonal entries of A_prime non-zero.
        # We must explicitly zero them out to match the definition of a simple graph.
        np.fill_diagonal(A_prime, 0)
        
        # Check if the transformation results in the same graph, i.e., T(G) = G.
        # This is true if their adjacency matrices are identical.
        if np.array_equal(A, A_prime):
            counts_by_n[n] = counts_by_n.get(n, 0) + 1

    # Format the output as an equation, as requested.
    if not counts_by_n:
        print("0 = 0")
        return

    # Sort the results by the number of vertices.
    sorted_counts = sorted(counts_by_n.items())
    
    # Create the equation string from the counts found.
    count_values = [str(count) for n, count in sorted_counts]
    total_count = sum(count for n, count in sorted_counts)
    equation_str = " + ".join(count_values)
    
    # Print the final equation.
    print(f"{equation_str} = {total_count}")

solve_graph_transformation()