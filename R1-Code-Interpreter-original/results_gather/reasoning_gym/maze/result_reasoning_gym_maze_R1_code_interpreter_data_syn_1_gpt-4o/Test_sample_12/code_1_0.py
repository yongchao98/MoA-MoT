from collections import deque

# Define the grid
grid = [
    "QQQQQQQ",
    "QQ%%%%Q",
    "QQ`%Q%Q",
    "Q%%Q%%Q",
    "Q%%%Q%Q",
    "Q%QQ%(Q",
    "QQQQQQQ"
]

# Define start and goal positions
start = (5, 5)
goal = (2, 2)

# Define possible moves: up, down, left, right
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS setup
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# BFS loop
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] == '%':
                queue.append(((nx, ny), steps + 1))
                visited.add((nx, ny))