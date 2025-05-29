from collections import deque

# Define the grid
grid = [
    ['-', '-', '-', 'p', 'q', 'p', 'p', '-'],
    ['-', '-', '-', '-', 'p', '-'],
    ['-', '-', '?', 'p', 'p', '-'],
    ['-', 'p', 'p', 'p', 'p', '-'],
    ['-', '-', '-', '-', '-', '-']
]

# Define the start and goal positions
start = (0, 3)
goal = (2, 2)

# Define possible movements: up, down, left, right
movements = [(-1, 0), (1, 0), (0, -1), (0, 1)]

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
        for dx, dy in movements:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] in ['p', '?']:  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)