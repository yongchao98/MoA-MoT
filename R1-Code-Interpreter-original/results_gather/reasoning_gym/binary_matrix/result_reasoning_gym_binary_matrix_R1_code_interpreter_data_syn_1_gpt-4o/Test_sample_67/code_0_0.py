from collections import deque

def nearest_zero_distance(matrix):
    n = len(matrix)
    dist = [[float('inf')] * n for _ in range(n)]
    queue = deque()

    # Initialize the queue with all 0 positions
    for i in range(n):
        for j in range(n):
            if matrix[i][j] == 0:
                queue.append((i, j))
                dist[i][j] = 0

    # Directions for moving in the grid (right, left, down, up)
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    # BFS
    while queue:
        x, y = queue.popleft()
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < n and 0 <= ny < n:
                if dist[nx][ny] > dist[x][y] + 1:
                    dist[nx][ny] = dist[x][y] + 1
                    queue.append((nx, ny))

    return dist

# Input matrix
matrix = [
    [1, 1, 1, 1, 0, 0, 0, 0, 0],
    [0, 1, 1, 1, 1, 0, 0, 1, 0],
    [0, 1, 1, 1, 1, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 1, 0, 1, 1],
    [1, 1, 1, 0, 1, 1, 1, 0, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 0],
    [1, 1, 1, 1, 0, 0, 0, 1, 0],
    [1, 1, 1, 0, 1, 1, 0, 0, 1],
    [1, 0, 1, 0, 1, 0, 1, 0, 1]
]

# Calculate the distance matrix
distance_matrix = nearest_zero_distance(matrix)

# Print the result
for row in distance_matrix:
    print(row)