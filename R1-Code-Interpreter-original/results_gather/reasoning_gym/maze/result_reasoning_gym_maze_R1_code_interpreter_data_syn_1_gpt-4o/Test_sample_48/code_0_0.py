from collections import deque

# Define the grid
grid = [
    "yyyyyy",
    "yyyy(y",
    "y{yOOy",
    "yOOOOy",
    "yOOOOy",
    "yyyyyy"
]

# Define the start and goal characters
start_char = '('
goal_char = '{'

# Find the start and goal positions
start_pos = None
goal_pos = None
for i, row in enumerate(grid):
    for j, char in enumerate(row):
        if char == start_char:
            start_pos = (i, j)
        elif char == goal_char:
            goal_pos = (i, j)

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

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
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] == 'O' or grid[nx][ny] == goal_char:
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Get the minimum number of steps
min_steps = bfs(start_pos, goal_pos)
print(min_steps)