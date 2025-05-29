from collections import deque

# Define the grid
grid = [
    "JJJJJJJJ",
    "J,,,J,,J",
    "JJJ,,,JJ",
    "J,J,,RJJ",
    "J,J,,J,J",
    "J,5,,,,J",
    "J,,,JJ,J",
    "JJJJJJJJ"
]

# Directions for moving up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the start position
start = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == '5':
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
    
    # Check if we reached the goal
    if grid[x][y] == 'R':
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] != 'J':  # Not a wall
                visited.add((nx, ny))
                queue.append((nx, ny, steps + 1))