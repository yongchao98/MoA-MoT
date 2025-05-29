from collections import deque

def find_min_steps(grid):
    # Locate start and goal positions
    start = None
    goal = None
    for r in range(len(grid)):
        for c in range(len(grid[0])):
            if grid[r][c] == ';':
                start = (r, c)
            elif grid[r][c] == 'J':
                goal = (r, c)
    
    # Directions for moving in the grid (up, down, left, right)
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    # BFS initialization
    queue = deque([(start[0], start[1], 0)])  # (row, col, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        r, c, steps = queue.popleft()
        
        # Check if we reached the goal
        if (r, c) == goal:
            return steps
        
        # Explore neighbors
        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if (0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and
                grid[nr][nc] in [']', 'J'] and (nr, nc) not in visited):
                visited.add((nr, nc))
                queue.append((nr, nc, steps + 1))
    
    return -1  # If no path is found

# Define the grid
grid = [
    "iiiiiiii",
    "i]ii]]ii",
    "i]iii]ii",
    "ii]]]J]i",
    "i]i]]]]i",
    "i]]]]]ii",
    "i;]i]i]i",
    "iiiiiiii"
]

# Find the minimum steps
min_steps = find_min_steps(grid)
print(min_steps)