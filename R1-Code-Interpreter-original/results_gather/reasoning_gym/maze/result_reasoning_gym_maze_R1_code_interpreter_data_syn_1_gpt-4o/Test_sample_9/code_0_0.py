from collections import deque

# Define the grid
grid = [
    "MMMMMMMM",
    "MGGtGMlM",
    "MMGMGGGM",
    "MMGGGGMM",
    "MGMMGMMM",
    "MGGGMGMM",
    "MGGGGGMM",
    "MMMMMMMM"
]

# Define the start and goal positions
start = (1, 2)  # 't' is at (1, 2)
goal = (1, 6)   # 'l' is at (1, 6)

# Define the possible movements (up, down, left, right)
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
                if grid[nx][ny] == 'G' or grid[nx][ny] == 'l':  # Can move to 'G' or 'l'
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)