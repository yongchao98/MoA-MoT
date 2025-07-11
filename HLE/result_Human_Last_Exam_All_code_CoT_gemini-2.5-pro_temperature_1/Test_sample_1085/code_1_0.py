import numpy as np

def solve_prime_path():
    """
    Calculates the number of distinct Prime Paths of length 4 from (1,1) to (5,7).
    """
    # The set of allowed coordinates based on the problem (1 and primes up to 7)
    coords = [1, 2, 3, 5, 7]
    n = len(coords)

    # Create a mapping from node coordinates to matrix indices
    node_map = {val: i for i, val in enumerate(coords)}
    
    # Calculate shortest path distance for explanation
    # The steps are the differences in indices in the coords list
    x_dist = node_map[5] - node_map[1]
    y_dist = node_map[7] - node_map[1]
    shortest_path_len = x_dist + y_dist
    required_path_len = 4

    # The number of paths is 0 if the required length is less than the shortest possible length.
    # We can stop here, but the code below will calculate it formally.
    if required_path_len < shortest_path_len:
        num_paths = 0
    else:
        # Formally calculate using graph theory for completeness
        # Total number of nodes in our relevant grid portion
        num_nodes = n * n 
        
        # Adjacency matrix: A[i,j] = 1 if a path exists between node i and j
        adj_matrix = np.zeros((num_nodes, num_nodes), dtype=int)
        
        for r_idx, r_val in enumerate(coords):
            for c_idx, c_val in enumerate(coords):
                node_idx = r_idx * n + c_idx
                
                # Horizontal adjacency (left/right)
                if c_idx > 0:
                    neighbor_idx = r_idx * n + (c_idx - 1)
                    adj_matrix[node_idx, neighbor_idx] = 1
                if c_idx < n - 1:
                    neighbor_idx = r_idx * n + (c_idx + 1)
                    adj_matrix[node_idx, neighbor_idx] = 1

                # Vertical adjacency (up/down)
                if r_idx > 0:
                    neighbor_idx = (r_idx - 1) * n + c_idx
                    adj_matrix[node_idx, neighbor_idx] = 1
                if r_idx < n - 1:
                    neighbor_idx = (r_idx + 1) * n + c_idx
                    adj_matrix[node_idx, neighbor_idx] = 1
        
        # The number of paths of length k is given by A^k
        # We use matrix_power to calculate A^4
        adj_matrix_pow_4 = np.linalg.matrix_power(adj_matrix, required_path_len)
        
        # Find the start and end indices in the matrix
        start_node_coord = (1, 1)
        end_node_coord = (5, 7)
        start_idx = node_map[start_node_coord[0]] * n + node_map[start_node_coord[1]]
        end_idx = node_map[end_node_coord[0]] * n + node_map[end_node_coord[1]]
        
        num_paths = adj_matrix_pow_4[start_idx, end_idx]

    # Final output as requested
    print("Final Equation:")
    print(f"Minimum horizontal moves (1->2->3->5) = {x_dist}")
    print(f"Minimum vertical moves (1->2->3->5->7) = {y_dist}")
    print(f"Shortest path length = {x_dist} + {y_dist} = {shortest_path_len}")
    print(f"Required path length = {required_path_len}")
    print(f"Since {required_path_len} < {shortest_path_len}, the number of paths is {num_paths}.")

solve_prime_path()