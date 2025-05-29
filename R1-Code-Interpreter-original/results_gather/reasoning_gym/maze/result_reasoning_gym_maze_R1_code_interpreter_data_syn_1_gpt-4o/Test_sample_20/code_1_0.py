from collections import deque

# Define the grid
grid = [
    "yyyyyyyy",
    "yy666yyy",
    "yyy666}y",
    "y666666y",
    "y666666y",
    "yy6666yy",
    "y66y6$yy",
    "yyyyyyyy"
]

# Find the start and goal positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == '}':
            # Start from the adjacent '6' cell
            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                ni, nj = i + dx, j + dy
                if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and grid[ni][nj] == '6':
                    start = (ni, nj)
                    break
        elif cell == '$':
            goal = (i, j)

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS to find the shortest path
def bfs(start, goal):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # Check if we have reached the goal
        if (x, y) == goal:
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] == '6' and (nx, ny) not in visited:
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Calculate the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)