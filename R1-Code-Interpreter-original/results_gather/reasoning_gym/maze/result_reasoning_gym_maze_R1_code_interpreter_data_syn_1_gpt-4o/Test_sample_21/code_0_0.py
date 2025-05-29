from collections import deque

# Define the grid
grid = [
    ";;;;;;;;;",
    ";;rrrr;;;",
    ";rrrrrrW;",
    ";;rrrrrr;",
    ";rrrrr;;;",
    ";rr;#rrr;",
    ";rr;rr;;;",
    ";r;rr;rr;",
    ";;;;;;;;;"
]

# Find the starting point 'W'
start = None
for i, row in enumerate(grid):
    if 'W' in row:
        start = (i, row.index('W'))
        break

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS to find the shortest path
def bfs(start):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # Check if we reached the goal
        if grid[x][y] == '#':
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] != ';' and (nx, ny) not in visited:
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))

# Get the minimum number of steps
min_steps = bfs(start)
print(min_steps)