from collections import deque

# Define the grid
grid = [
    ".........",
    ".IIIIIII.",
    ".I.IdI.I.",
    "...I.III.",
    ".II.III..",
    "..-III.I.",
    "..II.III.",
    ".IIIIII..",
    "........."
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Find the start and goal positions
start = None
goal = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == '-':
            start = (i, j)
        elif grid[i][j] == 'd':
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
        
        # Check if we reached the goal
        if (x, y) == goal:
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] == 'I' or grid[nx][ny] == 'd':  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Get the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)