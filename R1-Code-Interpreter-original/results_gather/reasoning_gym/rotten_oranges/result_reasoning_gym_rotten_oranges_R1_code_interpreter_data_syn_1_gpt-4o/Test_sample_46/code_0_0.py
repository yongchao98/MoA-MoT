from collections import deque

def orangesRotting(grid):
    rows, cols = len(grid), len(grid[0])
    queue = deque()
    fresh_oranges = 0

    # Initialize the queue with all rotten oranges and count fresh oranges
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 2:
                queue.append((r, c, 0))  # (row, col, time)
            elif grid[r][c] == 1:
                fresh_oranges += 1

    # Directions for 4-directional adjacency
    directions = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    minutes = 0

    # BFS
    while queue:
        r, c, time = queue.popleft()
        minutes = max(minutes, time)

        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols and grid[nr][nc] == 1:
                grid[nr][nc] = 2  # Make the orange rotten
                fresh_oranges -= 1
                queue.append((nr, nc, time + 1))

    # If there are still fresh oranges left, return -1
    if fresh_oranges > 0:
        return -1
    return minutes

# Define the grid
grid = [
    [1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0],
    [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
    [1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 0],
    [2, 1, 2, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1],
    [1, 2, 1, 1, 0, 1, 1, 2, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 1, 1, 0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1],
    [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1],
    [0, 1, 0, 1, 1, 0, 2, 0, 0, 1, 1, 2, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 2, 1, 1, 1],
    [1, 1, 1, 1, 2, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0],
    [1, 1, 1, 1, 0, 1, 1, 2, 2, 2, 1, 1, 1, 0, 1, 2],
    [0, 1, 1, 0, 1, 2, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0],
    [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 2],
    [1, 1, 0, 1, 1, 1, 0, 2, 1, 1, 1, 1, 0, 1, 1, 1],
    [1, 1, 1, 2, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1]
]

# Get the result
result = orangesRotting(grid)
print(result)