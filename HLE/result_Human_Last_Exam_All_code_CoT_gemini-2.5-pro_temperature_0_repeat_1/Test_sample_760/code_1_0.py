import networkx as nx
import numpy as np

def solve():
    """
    Finds all non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the transformation T(G) = G.
    """
    count = 0
    solutions_by_n = {i: 0 for i in range(1, 8)}
    found_graphs_info = []

    # nx.graph_atlas_g() returns a list of all non-isomorphic graphs up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        # We only consider connected graphs as per the problem statement.
        if not nx.is_connected(G):
            continue

        n = G.number_of_nodes()
        if n == 0:
            continue

        # Handle the trivial case of a single-vertex graph.
        # A=[[0]], A^2=[[0]]. T(G) has no edges. So T(G)=G.
        if n == 1:
            is_solution = True
        else:
            # Get the adjacency matrix A.
            A = nx.to_numpy_array(G, dtype=int)

            # Compute A^2. The entry (i, j) of A^2 is the number of paths of length 2.
            A2 = np.dot(A, A)

            # Construct the adjacency matrix A_prime for the transformed graph T(G).
            A_prime = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i + 1, n):
                    if 1 <= A2[i, j] <= 2:
                        A_prime[i, j] = 1
                        A_prime[j, i] = 1
            
            # Check if the transformation is an identity, i.e., T(G) = G.
            is_solution = np.array_equal(A, A_prime)

        if is_solution:
            count += 1
            solutions_by_n[n] += 1
            
            # Store info for printing later
            e = G.number_of_edges()
            name = "Unknown"
            if n == 1: name = "K1 (single vertex)"
            elif n == 3 and e == 3: name = "K3 (Triangle)"
            elif n == 4 and e == 6: name = "K4 (Complete graph)"
            elif n == 6 and e == 12: name = "Octahedral graph"
            
            found_graphs_info.append(f"- Graph with {n} vertices and {e} edges ({name})")

    # Print the results
    print("Found the following graphs that satisfy the condition T(G) = G:")
    for info in found_graphs_info:
        print(info)

    print("\nSummary by number of vertices:")
    equation_parts = []
    for n_val, num in solutions_by_n.items():
        if num > 0:
            print(f"  - {num} graph(s) with {n_val} vertices")
            equation_parts.append(str(num))

    # The prompt asks to output each number in the final equation.
    print("\nThe final count is the sum of these numbers:")
    equation = " + ".join(equation_parts)
    print(f"{equation} = {count}")

solve()
<<<4>>>