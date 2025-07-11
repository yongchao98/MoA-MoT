import networkx as nx
import numpy as np

def solve_graph_problem():
    """
    Finds and counts the non-isomorphic, connected graphs with at most 7 vertices
    that are fixed points of the specified transformation T.
    """
    # A dictionary to store the count of solutions for each number of vertices n.
    counts_by_n = {n: 0 for n in range(1, 8)}

    # nx.graph_atlas_g() provides an iterator over all non-isomorphic graphs
    # with up to 7 vertices. This is our source of candidate graphs.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        # The atlas yields an empty graph for n=0; we skip it.
        if n == 0:
            continue
        
        # Per the problem, we only consider connected graphs.
        if not nx.is_connected(G):
            continue
        
        # Obtain the adjacency matrix 'A' for the graph G.
        # We sort the nodes to ensure a consistent ordering for the matrix.
        A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))

        # Calculate A^2. The entry (A^2)_ij gives the number of length-2
        # paths between vertices i and j, which is also their number of common neighbors.
        A2 = np.dot(A, A)
        
        # Construct the adjacency matrix 'A_T' for the transformed graph T(G).
        # An edge (i, j) exists in T(G) if (A^2)_ij is 1 or 2.
        # We use numpy's `where` function for an efficient implementation of this rule.
        A_T = np.where((A2 == 1) | (A2 == 2), 1, 0)
        
        # The transformation does not create self-loops, so ensure the diagonal is zero.
        np.fill_diagonal(A_T, 0)
        
        # A graph G is a solution if T(G) = G, which means their adjacency matrices
        # must be identical.
        if np.array_equal(A, A_T):
            # If it is a solution, increment the counter for its number of vertices.
            counts_by_n[n] += 1

    # Prepare the final output as an equation showing the breakdown by n.
    count_values = list(counts_by_n.values())
    total_count = sum(count_values)
    equation_parts = [str(c) for c in count_values]
    
    # The numbers in the sum correspond to n=1, 2, 3, 4, 5, 6, and 7 respectively.
    equation = " + ".join(equation_parts) + f" = {total_count}"
    
    print(equation)

# Run the solver function
solve_graph_problem()