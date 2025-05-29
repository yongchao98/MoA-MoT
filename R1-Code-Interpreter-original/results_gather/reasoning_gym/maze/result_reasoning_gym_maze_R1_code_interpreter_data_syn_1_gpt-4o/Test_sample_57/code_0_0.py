from collections import deque

# Define the grid
grid = [
    "aaaaaaaaaa",
    "a:x:a:a::a",
    "a::::::{aa",
    "a::::::aaa",
    "a:a:a::::a",
    "a::::::a:a",
    "a::::a:a:a",
    "a::aaaa::a",
    "aa:aa::::a",
    "aaaaaaaaaa"
]

# Find the start and goal positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, char in enumerate(row):
        if char == '{':
            start = (i, j)
        elif char == 'x':
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
                if grid[nx][ny] == ':' or grid[nx][ny] == 'x':
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)