from collections import deque

# Define the grid
grid = [
    "qqqqqqq",
    "qZ||||q",
    "qq|||bq",
    "qq||q|q",
    "q|||||q",
    "qq||q|q",
    "qqqqqqq"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Find the start ('b') and goal ('Z') positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == 'b':
            start = (i, j)
        elif cell == 'Z':
            goal = (i, j)

# Directions for moving up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS to find the shortest path
def bfs(start, goal):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # If we reach the goal, return the number of steps
        if (x, y) == goal:
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] == '|':  # Can only move through passages
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Get the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)