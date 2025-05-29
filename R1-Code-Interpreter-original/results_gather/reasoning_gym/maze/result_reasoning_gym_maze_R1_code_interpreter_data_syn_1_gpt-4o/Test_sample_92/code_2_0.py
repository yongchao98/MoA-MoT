from collections import deque

# Define the grid
grid = [
    "zzzzzzzzzz",
    "z)z)))))zz",
    "z)m)))zz)z",
    "zzzz))))zz",
    "z)z)))zz)z",
    "zzzz)0zzzz",
    "z)z)z)z))z",
    "z)))z)))zz",
    "z)zz)))z)z",
    "zzzzzzzzzz"
]

# Correct start and goal positions
start = (5, 4)  # Position of '0'
goal = (2, 1)   # Position of 'm'

# Define possible movements: right, left, down, up
directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

# BFS to find the shortest path
def bfs(grid, start, goal):
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
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] == ')' or grid[nx][ny] == 'm':  # Check if it's a passage or the goal
                    queue.append(((nx, ny), steps + 1))
                    visited.add((nx, ny))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(grid, start, goal)
print(min_steps)