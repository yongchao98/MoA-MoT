from collections import deque

# Define the grid
grid = [
    "yyyyyyyyy",
    "yyRRRRyyy",
    "yRRyyRyyy",
    "yRRRRRRRy",
    "yR!yyRRyy",
    "yyRyyRRRy",
    "yyyyRRRyy",
    "yRRRwRRRy",
    "yyyyyyyyy"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Define directions for moving up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the start and goal positions
start = None
goal = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == '!':
            start = (i, j)
        elif grid[i][j] == 'w':
            goal = (i, j)

# BFS to find the shortest path
def bfs(start, goal):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # Check if we reached the goal
        if (x, y) == goal:
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] in ('R', 'w') and (nx, ny) not in visited:
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Get the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)