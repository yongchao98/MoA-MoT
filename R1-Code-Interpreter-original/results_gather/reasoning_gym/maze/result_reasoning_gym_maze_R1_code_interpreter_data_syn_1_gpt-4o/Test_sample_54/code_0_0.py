from collections import deque

# Define the grid
grid = [
    "333333333",
    "3****3**3",
    "333*333*3",
    "3**3*33*3",
    "313E****3",
    "3*3*3***3",
    "3****33*3",
    "3****3**3",
    "333333333"
]

# Find the start ('E') and goal ('1') positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == 'E':
            start = (i, j)
        elif cell == '1':
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
                grid[nx][ny] != '3' and (nx, ny) not in visited):
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))

# Calculate the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)