from collections import deque

# Define the grid
grid = [
    "&&&&&&&&&&",
    "&J3JJJ&&J&",
    "&&JJ&J&&J&",
    "&&&JJJ&-&&",
    "&&&&&JJJJ&",
    "&JJJJJ&J&&",
    "&&JJ&&JJ&&",
    "&JJJ&JJJ&&",
    "&&J&JJJ&&&",
    "&&&&&&&&&&"
]

# Define start and goal positions
start = (3, 7)
goal = (1, 1)

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS to find the shortest path
def bfs(start, goal, grid):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # Check if we reached the goal
        if (x, y) == goal:
            return steps
        
        # Explore neighbors
        for dx, dy in moves:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] != '&' and (nx, ny) not in visited:
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal, grid)
print(min_steps)