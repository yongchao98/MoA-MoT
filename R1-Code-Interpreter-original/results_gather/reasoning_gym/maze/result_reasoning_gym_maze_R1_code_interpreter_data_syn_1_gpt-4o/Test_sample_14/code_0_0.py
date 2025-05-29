from collections import deque

# Define the grid
grid = [
    "999999999",
    "999PP9PP9",
    "9999PPP99",
    "9P9PPP999",
    "9PPPP9P59",
    "9PPPPPPP9",
    "9PPP.9P99",
    "9P9PPPP99",
    "999999999"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Find the start position (5)
start_pos = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == '5':
            start_pos = (i, j)
            break
    if start_pos:
        break

# Define the directions for movement (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS initialization
queue = deque([(start_pos[0], start_pos[1], 0)])  # (row, col, steps)
visited = set()
visited.add(start_pos)

# BFS loop
while queue:
    x, y, steps = queue.popleft()
    
    # Check if we reached the goal
    if grid[x][y] == '.':
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] in ('P', '.'):
                visited.add((nx, ny))
                queue.append((nx, ny, steps + 1))