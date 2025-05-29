from collections import deque

# Define the grid
grid = [
    '"""""""""',
    '"n""nnnn"',
    '"nnnn*nn"',
    '"rnnnnnn"',
    '""nnn"n""',
    '"""nn"""',
    '"nn"nn"n"',
    '"nnnnnnn"',
    '"""""""""'
]

# Find the start and goal positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, char in enumerate(row):
        if char == '*':
            start = (i, j)
        elif char == 'r':
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
        
        # If we reach the goal, return the number of steps
        if (x, y) == goal:
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] == 'n' or grid[nx][ny] == 'r':
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Calculate the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)