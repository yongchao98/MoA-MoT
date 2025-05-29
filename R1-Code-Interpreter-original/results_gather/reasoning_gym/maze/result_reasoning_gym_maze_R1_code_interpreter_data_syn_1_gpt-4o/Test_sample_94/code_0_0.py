from collections import deque

# Define the grid
grid = [
    [',', ',', ',', ',', ',', ',', ',', ','],
    [',', '7', ',', '7', ',', '7', ',', ','],
    [',', '7', '7', '7', ',', ',', ',', ','],
    [',', ',', '7', '7', '7', '7', '/', ','],
    [',', 'W', '7', '7', ',', ',', '7', ','],
    [',', ',', '7', '7', '7', '7', ',', ','],
    [',', '7', ',', ',', ',', ',', ',', ','],
    [',', ',', ',', ',', ',', ',', ',', ',']
]

# Define the start and goal positions
start = (3, 6)
goal = (4, 1)

# Define the possible movements (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

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
                if grid[nx][ny] == '7' or grid[nx][ny] == 'W':  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(grid, start, goal)
print(min_steps)