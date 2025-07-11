import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy T(G) = G, where T is the specified transformation.
    """
    # Dictionary to store the count of graphs for each number of vertices
    counts_per_n = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() returns a list of all graphs with up to 7 vertices.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        # We are interested in graphs with at least 1 vertex
        if n == 0:
            continue

        # The graph must be connected
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A of the graph G
        # Using a sorted nodelist ensures a consistent matrix representation.
        A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))

        # Compute A^2. The entry (i,j) of A^2 gives the number of length-2 paths between i and j.
        A_squared = np.dot(A, A)

        # Construct the adjacency matrix A_prime for the transformed graph T(G)
        A_prime = np.zeros_like(A)
        for i in range(n):
            for j in range(i + 1, n):  # Iterate over pairs of distinct vertices
                num_paths = A_squared[i, j]
                # An edge exists in T(G) iff there are 1 or 2 paths of length 2
                if 1 <= num_paths <= 2:
                    A_prime[i, j] = 1
                    A_prime[j, i] = 1

        # Check if the transformation results in the same graph, i.e., T(G) = G
        if np.array_equal(A, A_prime):
            counts_per_n[n] += 1

    # Format the final output as an equation showing the breakdown by vertex count
    equation_parts = [str(counts_per_n[n]) for n in range(1, 8)]
    total_count = sum(counts_per_n.values())
    
    final_equation = " + ".join(equation_parts) + f" = {total_count}"
    
    print("The number of graphs G such that T(G)=G for n=1 to 7 is given by the sum:")
    print(final_equation)

solve()
<<<5>>>