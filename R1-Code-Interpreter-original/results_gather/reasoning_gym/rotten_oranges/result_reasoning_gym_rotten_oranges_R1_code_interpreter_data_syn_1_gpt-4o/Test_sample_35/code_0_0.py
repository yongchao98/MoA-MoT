from collections import deque

def orangesRotting(grid):
    rows, cols = len(grid), len(grid[0])
    queue = deque()
    fresh_oranges = 0

    # Initialize the queue with all rotten oranges and count fresh oranges
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 2:
                queue.append((r, c))
            elif grid[r][c] == 1:
                fresh_oranges += 1

    # Directions for adjacent cells (up, down, left, right)
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    minutes_passed = 0

    # BFS to rot adjacent fresh oranges
    while queue and fresh_oranges > 0:
        minutes_passed += 1
        for _ in range(len(queue)):
            x, y = queue.popleft()
            for dx, dy in directions:
                nx, ny = x + dx, y + dy
                if 0 <= nx < rows and 0 <= ny < cols and grid[nx][ny] == 1:
                    grid[nx][ny] = 2
                    fresh_oranges -= 1
                    queue.append((nx, ny))

    # If there are still fresh oranges left, return -1
    return minutes_passed if fresh_oranges == 0 else -1

# Define the grid
grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1],
    [0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1],
    [1, 1, 1, 2, 1, 1, 1, 1, 2, 0, 1, 2],
    [1, 1, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1],
    [0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1],
    [1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 2, 1],
    [2, 2, 1, 1, 1, 1, 0, 1, 1, 1, 2, 2],
    [1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 2, 1]
]

# Calculate the minimum number of minutes
result = orangesRotting(grid)
print(result)