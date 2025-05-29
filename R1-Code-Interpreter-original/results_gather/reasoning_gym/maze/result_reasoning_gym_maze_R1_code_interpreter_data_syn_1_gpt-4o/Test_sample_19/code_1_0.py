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

# Corrected start and goal positions
start = (3, 5)  # Corrected position of ')'
goal = (6, 5)   # Position of 'k'

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
            if 0 <= nx < 9 and 0 <= ny < 9 and grid[nx][ny] == '#' and (nx, ny) not in visited:
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)