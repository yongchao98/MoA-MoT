from collections import deque

def min_steps_to_goal(grid):
    # Directions for moving up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    # Find the start position
    start = None
    for i, row in enumerate(grid):
        for j, value in enumerate(row):
            if value == '#':
                start = (i, j)
                break
        if start:
            break
    
    # Initialize BFS
    queue = deque([(start[0], start[1], 0)])  # (row, col, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        x, y, steps = queue.popleft()
        
        # Check if we reached the goal
        if grid[x][y] == 'Z':
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] in ('f', 'Z'):  # Can move to passage or goal
                    queue.append((nx, ny, steps + 1))
                    visited.add((nx, ny))
    
    return -1  # If no path is found

# Define the grid
grid = [
    "zzzzzz",
    "zffffz",
    "z#fzZz",
    "zfffzz",
    "zffzzz",
    "zzzzzz"
]

# Calculate the minimum steps
min_steps = min_steps_to_goal(grid)
print(min_steps)