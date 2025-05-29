from collections import deque

# Define the grid
grid = [
    "--------",
    "-ooo-o--",
    "-Jooooo-",
    "-o-o-(o-",
    "-ooooo--",
    "-o-oo-o-",
    "--oo-oo-",
    "--------"
]

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the start and goal positions
start = None
goal = None
for r in range(len(grid)):
    for c in range(len(grid[0])):
        if grid[r][c] == '(':
            start = (r, c)
        elif grid[r][c] == 'J':
            goal = (r, c)

# BFS to find the shortest path
def bfs(start, goal):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (current, steps) = queue.popleft()
        
        if current == goal:
            return steps
        
        for direction in directions:
            nr, nc = current[0] + direction[0], current[1] + direction[1]
            if (0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and
                (grid[nr][nc] == 'o' or (nr, nc) == goal) and (nr, nc) not in visited):
                visited.add((nr, nc))
                queue.append(((nr, nc), steps + 1))
    
    return -1  # If no path is found

# Calculate the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)