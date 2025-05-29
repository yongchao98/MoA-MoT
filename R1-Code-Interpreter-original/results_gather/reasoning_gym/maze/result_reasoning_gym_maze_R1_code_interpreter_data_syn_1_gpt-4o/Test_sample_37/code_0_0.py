from collections import deque

# Define the grid
grid = [
    "^^^^^^^^",
    "^4^44c4^",
    "^44^4^4^",
    "^4^^4^4^",
    "^44l4^^^",
    "^4^44^^^",
    "^4^^44^^",
    "^^^^^^^^"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Directions for moving in the grid (right, down, left, up)
directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]

# Find the start position 'c'
start = None
for i in range(len(grid)):
    for j in range(len(grid[0])):
        if grid[i][j] == 'c':
            start = (i, j)
            break
    if start:
        break

# BFS setup
queue = deque([(start[0], start[1], 0)])  # (row, col, steps)
visited = set()
visited.add(start)

# BFS loop
while queue:
    x, y, steps = queue.popleft()
    
    # Check if we've reached the goal 'l'
    if grid[x][y] == 'l':
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] != '^':  # Check if it's not a wall
                visited.add((nx, ny))
                queue.append((nx, ny, steps + 1))