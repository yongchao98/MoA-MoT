from collections import deque

# Define the grid
grid = [
    ">>>>>>>>",
    ">RRRRS>>",
    ">>>RRRRR",
    ">RRRRR>R",
    ">RRRRRR}>",
    ">RR>R>>R>",
    ">RR>>>>>",
    ">>RRRR>R>",
    ">>>>>>>>"
]

# Find the start position 'S'
start = None
for i, row in enumerate(grid):
    if 'S' in row:
        start = (i, row.index('S'))
        break

# Define the directions for movement (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS to find the shortest path
def bfs(start):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # Check if we have reached the goal
        if grid[x][y] == '}':
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] == 'R' or grid[nx][ny] == '}':
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start)
print(min_steps)