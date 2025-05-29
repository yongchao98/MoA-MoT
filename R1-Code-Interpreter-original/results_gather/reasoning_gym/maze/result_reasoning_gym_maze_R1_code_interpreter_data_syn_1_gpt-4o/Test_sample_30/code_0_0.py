from collections import deque

# Define the grid
grid = [
    [',', ',', ',', ',', ',', ',', ',', ','],
    [',', 'm', 'm', 'm', 'm', 'm', 'm', ','],
    [',', 'm', ',', 'm', 'm', 'm', ',', ','],
    [',', ',', 'u', 'm', ',', 'm', 'm', ','],
    [',', 'm', 'm', ',', ',', 'm', ',', ','],
    [',', 'm', ',', 'm', 'm', 'm', 'm', ','],
    [',', 'm', 'm', '8', 'm', 'm', 'm', ','],
    [',', ',', ',', ',', ',', ',', ',', ',']
]

# Locate start and goal
start = (6, 3)
goal = (3, 2)

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
                if grid[nx][ny] == 'm' or grid[nx][ny] == 'u':  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)