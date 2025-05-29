from collections import deque

# Define the grid
grid = [
    "wwwwwww",
    "wwpwwww",
    "wwFFFww",
    "wwFwFuw",
    "wFFwFww",
    "wwFFFFw",
    "wwwwwww"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Define directions for moving up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the starting point 'p'
start = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == 'p':
            start = (i, j)
            break
    if start:
        break

# Initialize the queue for BFS
queue = deque([(start[0], start[1], 0)])  # (row, col, steps)
visited = set()
visited.add(start)

# Perform BFS
while queue:
    x, y, steps = queue.popleft()
    
    # Check if we have reached the goal 'u'
    if grid[x][y] == 'u':
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] in ('F', 'u'):  # Can move to passage or goal
                queue.append((nx, ny, steps + 1))
                visited.add((nx, ny))