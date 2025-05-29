from collections import deque

# Define the grid
grid = [
    "{{{{{{{{{",
    "{{Z{{Z{Z{",
    "{ZZZ{Z2Z{",
    "{ZZZZ{ZZ{",
    "{Z{ZZ{ZZ{",
    "{ZZ*ZZZZ{",
    "{Z{ZZZZ{{",
    "{{{{ZZZZ{",
    "{{{{{{{{{"
]

# Define the start and goal characters
start_char = '2'
goal_char = '*'

# Find the start position
start_pos = None
for i, row in enumerate(grid):
    if start_char in row:
        start_pos = (i, row.index(start_char))
        break

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS to find the shortest path
def bfs(start, grid):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # Check if we reached the goal
        if grid[x][y] == goal_char:
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] in ('Z', goal_char):
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If goal is not reachable

# Find the minimum number of steps
min_steps = bfs(start_pos, grid)
print(min_steps)