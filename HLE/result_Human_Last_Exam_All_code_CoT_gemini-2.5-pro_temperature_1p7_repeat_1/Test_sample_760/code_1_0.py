import networkx as nx
import numpy as np

def solve():
    """
    This function finds the number of non-isomorphic, connected graphs G 
    with at most 7 vertices that are fixed points of the transformation T,
    i.e., T(G) = G.
    """
    
    # counts[n] will store the number of solution graphs for n vertices.
    # Index 0 is unused; indices 1-7 correspond to graphs with 1-7 vertices.
    counts = [0] * 8

    # networkx.graph_atlas_g() provides a list of all non-isomorphic graphs
    # with up to 7 vertices.
    graph_atlas = nx.graph_atlas_g()

    for G in graph_atlas:
        n = G.number_of_nodes()

        # We only consider connected graphs as per the problem statement.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A for the graph.
        # We sort the nodelist to ensure a canonical matrix representation.
        nodes = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodes, dtype=np.int32)

        # For a graph with 1 vertex, T(G)=G. It is connected, so it is a solution.
        if n == 1:
            counts[1] += 1
            continue
        
        # For n=0, is_connected is false, so it's already skipped.
        # For n > 1, we perform the check using matrix operations.
        if n > 1:
            # Compute A^2. The entry (i, j) of A^2 gives the number of length-2 paths.
            A2 = A @ A

            # Construct the adjacency matrix A_prime for the transformed graph T(G).
            A_prime = np.zeros_like(A, dtype=np.int32)
            for i in range(n):
                for j in range(i + 1, n):
                    # An edge exists in T(G) if the number of length-2 paths is 1 or 2.
                    if A2[i, j] == 1 or A2[i, j] == 2:
                        A_prime[i, j] = 1
                        A_prime[j, i] = 1
            
            # Check if the transformation results in the same graph (A == A_prime).
            if np.array_equal(A, A_prime):
                counts[n] += 1

    total_count = sum(counts)
    
    # As requested, output the components of the sum to form an equation.
    equation_parts = []
    # Loop from n=1 to n=7
    for n in range(1, len(counts)):
        if counts[n] > 0:
            equation_parts.append(str(counts[n]))

    equation_str = " + ".join(equation_parts)
    
    print("The final count is the sum of the number of solution graphs found for each vertex count:")
    print(f"Number of solutions per vertex count (n=1, 3, 4, 6): {', '.join(equation_parts)}")
    print(f"Final Equation: {equation_str} = {total_count}")


solve()
<<<4>>>