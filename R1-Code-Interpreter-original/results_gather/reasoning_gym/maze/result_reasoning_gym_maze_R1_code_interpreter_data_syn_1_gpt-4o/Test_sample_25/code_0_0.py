from collections import deque

# Define the grid
grid = [
    "5555555555",
    "5LLLLLLL55",
    "555LLLLLL5",
    "55LLLLL.L5",
    "5LL55LLLL5",
    "5LL5LLLLL5",
    "55L5LLLL55",
    "5LLL:L5LL5",
    "5L55LLLLL5",
    "5555555555"
]

# Find the start ('.') and goal (':') positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == '.':
            start = (i, j)
        elif cell == ':':
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
            if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and
                grid[nx][ny] in 'L:' and (nx, ny) not in visited):
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Calculate the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)