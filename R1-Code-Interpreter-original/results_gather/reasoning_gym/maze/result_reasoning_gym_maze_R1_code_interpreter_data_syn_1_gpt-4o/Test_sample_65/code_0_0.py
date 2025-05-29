from collections import deque

# Define the grid
grid = [
    "ooooooooo",
    "oSooSSSoo",
    "ooSoooSoo",
    "oooSSSoSo",
    "ooSSooSoo",
    "oSSLoSSSo",
    "oSSSSSo6o",
    "oSSSoSoSo",
    "ooooooooo"
]

# Define directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the starting point 'L' and the goal '6'
start = None
goal = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == 'L':
            start = (i, j)
        elif cell == '6':
            goal = (i, j)

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
            if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and
                grid[nx][ny] == 'S' and (nx, ny) not in visited):
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)