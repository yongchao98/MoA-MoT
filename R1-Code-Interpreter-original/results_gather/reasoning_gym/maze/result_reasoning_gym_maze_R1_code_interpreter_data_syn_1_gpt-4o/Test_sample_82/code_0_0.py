from collections import deque

# Define the grid
grid = [
    "uuuuuuuu",
    "u+uu+++u",
    "u+++++uu",
    "uu+u+uuu",
    "u++u++Hu",
    "u+++uu+u",
    "u+c++++u",
    "uuuuuuuu"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Find the start position 'c'
start = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == 'c':
            start = (i, j)
            break
    if start:
        break

# BFS setup
queue = deque([(start[0], start[1], 0)])  # (row, col, steps)
visited = set()
visited.add(start)

# Direction vectors for moving up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Perform BFS
while queue:
    x, y, steps = queue.popleft()
    
    # Check if we've reached the goal
    if grid[x][y] == 'H':
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] in ('+', 'H'):  # Can move to passage or goal
                visited.add((nx, ny))
                queue.append((nx, ny, steps + 1))