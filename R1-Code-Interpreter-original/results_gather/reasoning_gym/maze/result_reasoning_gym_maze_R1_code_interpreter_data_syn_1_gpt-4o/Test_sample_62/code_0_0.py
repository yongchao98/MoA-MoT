from collections import deque

# Define the grid
grid = [
    "(((((((((",
    "(n(((nnn(",
    "(nn(nnn((",
    "(nnnn(nn(",
    "(f(nnn(((",
    "(nnn(nn((",
    "((nnnnnn(",
    "((nnnnnz(",
    "(((((((((("
]

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the start and goal positions
start = None
goal = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == 'f':
            start = (i, j)
        elif grid[i][j] == 'z':
            goal = (i, j)

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
                if grid[nx][ny] == 'n' or grid[nx][ny] == 'z':  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)