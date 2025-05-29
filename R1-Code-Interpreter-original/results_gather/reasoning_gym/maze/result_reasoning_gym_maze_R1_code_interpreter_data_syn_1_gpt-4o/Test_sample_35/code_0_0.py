from collections import deque

# Define the grid
grid = [
    "999999999",
    "9eee9e9e9",
    "99Xeeeee9",
    "9eee9ee99",
    "9e9ee9ee9",
    "9ee99eee9",
    "9eeee9ee9",
    "9eeQ9e999",
    "999999999"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the starting point 'Q'
start = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == 'Q':
            start = (i, j)
            break
    if start:
        break

# BFS initialization
queue = deque([(start[0], start[1], 0)])  # (row, col, steps)
visited = set()
visited.add(start)

# BFS loop
while queue:
    x, y, steps = queue.popleft()
    
    # Check if we reached the goal 'X'
    if grid[x][y] == 'X':
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] in ('e', 'X'):  # Can move to passage or goal
                visited.add((nx, ny))
                queue.append((nx, ny, steps + 1))