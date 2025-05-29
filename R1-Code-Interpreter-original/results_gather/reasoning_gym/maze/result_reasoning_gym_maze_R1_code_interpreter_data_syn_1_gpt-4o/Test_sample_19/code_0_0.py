from collections import deque

# Define the grid
grid = [
    "vvvvvvvvv",
    "vv#vvvvvv",
    "v###vvvvv",
    "v####)#vv",
    "v#v##v##v",
    "v##v#v##v",
    "v###vk#vv",
    "v###v##vv",
    "vvvvvvvvv"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Find the start position
start = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == ')':
            start = (i, j)
            break
    if start:
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
        
        # Check if we've reached the goal
        if grid[x][y] == 'k':
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] == '#' or grid[nx][ny] == 'k':
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start)
print(min_steps)