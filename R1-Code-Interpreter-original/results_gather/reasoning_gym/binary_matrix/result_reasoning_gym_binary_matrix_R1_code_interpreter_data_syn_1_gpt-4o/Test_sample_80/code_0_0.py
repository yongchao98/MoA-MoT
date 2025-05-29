from collections import deque

def nearest_zero_distance(matrix):
    rows, cols = len(matrix), len(matrix[0])
    distance = [[float('inf')] * cols for _ in range(rows)]
    queue = deque()

    # Initialize the queue with all 0 positions and set their distance to 0
    for r in range(rows):
        for c in range(cols):
            if matrix[r][c] == 0:
                queue.append((r, c))
                distance[r][c] = 0

    # Directions for moving up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    # BFS to calculate the minimum distance to a 0 for each cell
    while queue:
        r, c = queue.popleft()
        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                if distance[nr][nc] > distance[r][c] + 1:
                    distance[nr][nc] = distance[r][c] + 1
                    queue.append((nr, nc))

    return distance

# Input matrix
matrix = [
    [1, 1, 1, 1, 0, 0, 1, 1, 0, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 0, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 0, 1, 0, 1],
    [0, 0, 0, 1, 1, 1, 1, 0, 1, 1],
    [1, 1, 0, 1, 1, 0, 1, 0, 1, 1],
    [1, 1, 0, 1, 1, 1, 1, 0, 1, 1],
    [1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
    [1, 1, 1, 0, 1, 1, 0, 1, 1, 1]
]

# Calculate the distance matrix
distance_matrix = nearest_zero_distance(matrix)

# Print the result
print(distance_matrix)