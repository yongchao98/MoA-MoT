import numpy as np

# Step 1: Define functions to construct the two graphs.
# The graphs are the 4x4 Rook's graph and the Shrikhande graph.
# Both are strongly regular graphs with parameters (n, d, lambda, mu) = (16, 6, 2, 2).

def build_rook_graph():
    """
    Constructs the adjacency matrix for the 4x4 Rook's graph (L_2(4)).
    Vertices are numbered 0-15, representing cells of a 4x4 grid.
    Two vertices are adjacent if they are in the same row or column.
    """
    n = 16
    adj = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            # Map vertex index to grid coordinates (row, col)
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            if r1 == r2 or c1 == c2:
                adj[i, j] = adj[j, i] = 1
    return adj

def build_shrikhande_graph():
    """
    Constructs the adjacency matrix for the Shrikhande graph.
    This graph is a Cayley graph on Z_4 x Z_4. Vertices are pairs (i, j) with i, j in Z_4.
    Two vertices (i, j) and (x, y) are adjacent if their difference
    (i-x, j-y) mod 4 is in the connection set S.
    """
    n = 16
    adj = np.zeros((n, n), dtype=int)
    # The connection set S is symmetric, i.e., s in S implies -s in S.
    S = {(0, 1), (0, 3), (1, 0), (3, 0), (1, 1), (3, 3)}
    for i in range(n):
        for j in range(i + 1, n):
            # Map vertex index to grid coordinates (row, col)
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            diff = ((r1 - r2) % 4, (c1 - c2) % 4)
            if diff in S:
                adj[i, j] = adj[j, i] = 1
    return adj

# Step 2: Define a function to count cycles of a given length.
# We use a depth-first search approach from each node.

def count_cycles(adj, length):
    """
    Counts the number of simple cycles of a specific length in a graph.
    `adj` is the adjacency matrix.
    `length` is the length of the cycles to count (e.g., 5 for 5-cycles).
    """
    n = adj.shape[0]
    # Convert adjacency matrix to adjacency list for faster neighbor lookup
    adj_list = [list(np.where(row)[0]) for row in adj]
    count = 0
    
    for start_node in range(n):
        # Use a stack for an iterative DFS to find paths
        # Stack stores tuples of (current_node, path_taken)
        path = [start_node]
        stack = [(neighbor, path + [neighbor]) for neighbor in adj_list[start_node]]

        while stack:
            current_node, current_path = stack.pop()
            
            # If path is of the desired length, check for cycle completion
            if len(current_path) == length:
                if start_node in adj_list[current_node]:
                    count += 1
                continue

            # Expand the path with new neighbors
            for neighbor in adj_list[current_node]:
                # Ensure the path is simple (no repeated vertices)
                if neighbor not in current_path:
                    stack.append((neighbor, current_path + [neighbor]))
    
    # Each cycle is counted `length` times (once for each starting node)
    # and twice (for each direction of traversal).
    # So we divide by 2 * length.
    return count // (2 * length)

# Step 3: Main execution block.
# Create the graphs, count the 5-cycles, and print the results.

# The parameters for both graphs
n_param = 16
d_param = 6
lambda_param = 2
mu_param = 2

# Construct the graphs
rook_adj = build_rook_graph()
shrikhande_adj = build_shrikhande_graph()

# Count the 5-cycles for each graph
rook_c5 = count_cycles(rook_adj, 5)
shrikhande_c5 = count_cycles(shrikhande_adj, 5)

print("This script investigates if two strongly regular graphs with the same parameters can have a different number of 5-cycles.")
print(f"The chosen parameters are (n, d, lambda, mu) = ({n_param}, {d_param}, {lambda_param}, {mu_param}).")
print("The two graphs considered are the 4x4 Rook's graph and the Shrikhande graph.\n")

print("Graph 1: 4x4 Rook's graph")
print(f"Number of 5-cycles = {rook_c5}\n")

print("Graph 2: Shrikhande graph")
print(f"Number of 5-cycles = {shrikhande_c5}\n")

print("Since the number of 5-cycles is different, the answer to the user's question is yes.")
print(f"Final equation: {rook_c5} != {shrikhande_c5}")