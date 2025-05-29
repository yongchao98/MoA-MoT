from collections import deque

# Define the grid
grid = [
    "JJJJJJJJ",
    "J8888J8J",
    "JJ8JJJ8J",
    "JJ88888J",
    "J88J8)8J",
    "J8JJJ88J",
    "JT8J8J8J",
    "JJJJJJJJ"
]

# Define directions for moving up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the start and goal positions
start = None
goal = None
for r in range(len(grid)):
    for c in range(len(grid[0])):
        if grid[r][c] == 'T':
            start = (r, c)
        elif grid[r][c] == ')':
            goal = (r, c)

# BFS to find the shortest path
def bfs(start, goal):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (current_r, current_c), steps = queue.popleft()
        
        # Check if we have reached the goal
        if (current_r, current_c) == goal:
            return steps
        
        # Explore neighbors
        for dr, dc in directions:
            new_r, new_c = current_r + dr, current_c + dc
            if (0 <= new_r < len(grid) and 0 <= new_c < len(grid[0]) and
                grid[new_r][new_c] != 'J' and (new_r, new_c) not in visited):
                visited.add((new_r, new_c))
                queue.append(((new_r, new_c), steps + 1))
    
    return -1  # If no path is found

# Get the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)