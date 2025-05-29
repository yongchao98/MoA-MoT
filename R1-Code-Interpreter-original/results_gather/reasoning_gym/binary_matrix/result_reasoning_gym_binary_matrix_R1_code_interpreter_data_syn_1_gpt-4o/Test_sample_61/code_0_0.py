from collections import deque

def nearest_zero_distance(matrix):
    rows, cols = len(matrix), len(matrix[0])
    distances = [[float('inf')] * cols for _ in range(rows)]
    queue = deque()

    # Initialize the queue with all 0 positions
    for r in range(rows):
        for c in range(cols):
            if matrix[r][c] == 0:
                queue.append((r, c))
                distances[r][c] = 0

    # Directions for moving in the grid (up, down, left, right)
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    # Perform BFS
    while queue:
        r, c = queue.popleft()
        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols and distances[nr][nc] == float('inf'):
                distances[nr][nc] = distances[r][c] + 1
                queue.append((nr, nc))

    return distances

# Input matrix
matrix = [
    [0, 0, 1, 1, 0, 0, 0, 1, 0, 0],
    [1, 1, 1, 1, 1, 1, 0, 1, 1, 1],
    [0, 1, 1, 1, 0, 1, 0, 1, 0, 1],
    [1, 1, 0, 1, 1, 0, 1, 0, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 0, 1, 0, 0, 1, 1],
    [0, 1, 1, 1, 0, 1, 1, 1, 0, 1],
    [1, 1, 0, 1, 1, 1, 1, 1, 1, 0],
    [0, 1, 1, 1, 1, 0, 1, 1, 1, 1],
    [1, 0, 1, 0, 1, 0, 1, 1, 1, 1]
]

# Calculate distances
output_matrix = nearest_zero_distance(matrix)

# Print the output matrix
for row in output_matrix:
    print(row)