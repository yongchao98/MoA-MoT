from collections import deque

def nearest_zero_distance(matrix):
    n = len(matrix)
    m = len(matrix[0])
    distance = [[float('inf')] * m for _ in range(n)]
    queue = deque()

    # Initialize the queue with all 0 positions and set their distance to 0
    for i in range(n):
        for j in range(m):
            if matrix[i][j] == 0:
                queue.append((i, j))
                distance[i][j] = 0

    # Directions for moving up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    # BFS
    while queue:
        x, y = queue.popleft()
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < n and 0 <= ny < m:
                if distance[nx][ny] > distance[x][y] + 1:
                    distance[nx][ny] = distance[x][y] + 1
                    queue.append((nx, ny))

    return distance

# Input matrix
matrix = [
    [0, 1, 0, 1, 0, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 0, 1],
    [1, 1, 1, 1, 1, 1, 0, 1, 0],
    [0, 0, 1, 1, 1, 1, 0, 1, 1],
    [1, 1, 1, 0, 1, 1, 1, 0, 1],
    [1, 1, 1, 1, 0, 0, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 0, 1, 0],
    [0, 1, 0, 1, 0, 1, 1, 0, 0],
    [0, 1, 1, 1, 1, 0, 1, 0, 1]
]

# Calculate the distance matrix
result = nearest_zero_distance(matrix)

# Print the result
for row in result:
    print(row)