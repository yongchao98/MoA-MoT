from collections import deque

def min_steps_to_goal(grid):
    # Directions for moving up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    # Find the start and goal positions
    start = None
    goal = None
    for i, row in enumerate(grid):
        for j, cell in enumerate(row):
            if cell == 'G':
                start = (i, j)
            elif cell == 'n':
                goal = (i, j)
    
    # Initialize BFS
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
                if grid[nx][ny] in ('o', 'n'):  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Define the grid
grid = [
    "OOOOOOOOOO",
    "OoOooooOoO",
    "OooooooOoO",
    "OoOoooooOO",
    "OoooOOooGO",
    "OonoOOoOoO",
    "OooooooOOO",
    "OooooooOOO",
    "OooooooooO",
    "OOOOOOOOOO"
]

# Convert grid to list of lists
grid = [list(row) for row in grid]

# Get the minimum steps
min_steps = min_steps_to_goal(grid)
print(min_steps)