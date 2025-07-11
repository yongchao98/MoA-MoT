import numpy as np

def create_rooks_graph(k=4):
    """
    Creates the adjacency matrix for a k x k Rook's graph.
    Vertices are cells (i, j) in a k x k grid.
    Two vertices are adjacent if they are in the same row or column.
    """
    n = k * k
    adj_matrix = np.zeros((n, n), dtype=int)
    for v1 in range(n):
        i1, j1 = v1 // k, v1 % k
        for v2 in range(v1 + 1, n):
            i2, j2 = v2 // k, v2 % k
            if i1 == i2 or j1 == j2:
                adj_matrix[v1, v2] = 1
                adj_matrix[v2, v1] = 1
    return adj_matrix

def create_shrikhande_graph():
    """
    Creates the adjacency matrix for the Shrikhande graph.
    It's the Cayley graph of Z_4 x Z_4 with a specific connection set.
    """
    n = 16
    k = 4
    adj_matrix = np.zeros((n, n), dtype=int)
    # Connection set S = {(0,1), (0,3), (1,0), (3,0), (1,1), (3,3)}
    # or { (0,+-1), (+-1,0), (+-1,+-1) }
    connection_set = {(0, 1), (0, 3), (1, 0), (3, 0), (1, 1), (3, 3)}
    for v1 in range(n):
        i1, j1 = v1 // k, v1 % k
        for dr, dc in connection_set:
            i2 = (i1 + dr) % k
            j2 = (j1 + dc) % k
            v2 = i2 * k + j2
            adj_matrix[v1, v2] = 1
    return adj_matrix

def count_c5(A):
    """
    Counts the number of 5-cycles in a graph with adjacency matrix A.
    """
    n = A.shape[0]
    count = 0
    # This method finds paths v1-v2-v3-v4-v5 and checks if v5 is connected to v1.
    # It's designed to be clear and correct for small n.
    for v1 in range(n):
        # Find paths of length 2 starting at v1: v1-v2-v3
        neighbors_v1 = np.where(A[v1, :])[0]
        for v2 in neighbors_v1:
            neighbors_v2 = np.where(A[v2, :])[0]
            for v3 in neighbors_v2:
                if v3 == v1:
                    continue
                # Find paths of length 3 starting at v1: v1-v2-v3-v4
                neighbors_v3 = np.where(A[v3, :])[0]
                for v4 in neighbors_v3:
                    if v4 == v2 or v4 == v1:
                        continue
                    # Find paths of length 4 starting at v1: v1-v2-v3-v4-v5
                    neighbors_v4 = np.where(A[v4, :])[0]
                    for v5 in neighbors_v4:
                        if v5 == v3 or v5 == v2:
                            continue
                        # Check if v5 is connected to v1 to close the cycle,
                        # and ensure v5 is not v1 itself.
                        if A[v5, v1] and v5 != v1:
                            count += 1
    
    # Each 5-cycle is counted for each of its 5 vertices in each of 2 directions.
    # So, each cycle is counted 5 * 2 = 10 times.
    return count // 10

def main():
    """
    Main function to perform the calculation and print the result.
    """
    # Parameters for the strongly regular graphs
    n, d, lam, mu = 16, 6, 2, 2
    
    print(f"We investigate strongly regular graphs with parameters (n,d,lambda,mu) = ({n},{d},{lam},{mu}).")
    print("-" * 40)

    # Create the two graphs
    rooks_graph = create_rooks_graph(k=4)
    shrikhande_graph = create_shrikhande_graph()

    # Count the 5-cycles in each graph
    c5_rooks = count_c5(rooks_graph)
    c5_shrikhande = count_c5(shrikhande_graph)

    print(f"The 4x4 Rook's graph (L2(4)) has {c5_rooks} 5-cycles.")
    print(f"The Shrikhande graph has {c5_shrikhande} 5-cycles.")
    print("-" * 40)
    
    print("The final equation demonstrating they have a different number of 5-cycles is:")
    print(f"{c5_rooks} != {c5_shrikhande}")

if __name__ == "__main__":
    main()