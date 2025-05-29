from collections import deque

def find_minimum_steps(grid):
    # Directions for moving up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    # Find the start and goal positions
    start = None
    goal = None
    for i, row in enumerate(grid):
        for j, cell in enumerate(row):
            if cell == 'w':
                start = (i, j)
            elif cell == 'Z':
                goal = (i, j)
    
    # Initialize the queue for BFS
    queue = deque([(start[0], start[1], 0)])  # (row, col, steps)
    visited = set()
    visited.add(start)
    
    # Perform BFS
    while queue:
        x, y, steps = queue.popleft()
        
        # Check if we have reached the goal
        if (x, y) == goal:
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] == 'B' or grid[nx][ny] == 'Z':
                    visited.add((nx, ny))
                    queue.append((nx, ny, steps + 1))
    
    return -1  # If no path is found

# Define the grid
grid = [
    "UUUUUU",
    "UBUUUU",
    "UBBBBU",
    "UwUBUU",
    "UBUZBU",
    "UUUUUU"
]

# Find the minimum number of steps
min_steps = find_minimum_steps(grid)
print(min_steps)