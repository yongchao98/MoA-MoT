from collections import deque

# Define the grid
grid = [
    "qqqqqqqqq",
    "qqNxxxxxq",
    "qxxqxxxqq",
    "qxxqx3xqq",
    "qqxxqxxxq",
    "qxxxxxxxq",
    "qqqxxxxxq",
    "qxxqqxxxq",
    "qqqqqqqqq"
]

# Find the start position
start = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == 'N':
            start = (i, j)
            break
    if start:
        break

# Define the goal
goal = '3'

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS initialization
queue = deque([(start[0], start[1], 0)])  # (row, col, steps)
visited = set()
visited.add(start)

# BFS loop
while queue:
    x, y, steps = queue.popleft()
    
    # Check if we reached the goal
    if grid[x][y] == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] != 'q':  # Check if it's not a wall
                visited.add((nx, ny))
                queue.append((nx, ny, steps + 1))