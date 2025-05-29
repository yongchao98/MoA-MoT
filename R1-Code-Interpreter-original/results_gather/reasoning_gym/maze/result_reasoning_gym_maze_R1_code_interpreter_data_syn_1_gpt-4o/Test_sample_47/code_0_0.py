from collections import deque

def min_steps_to_goal(grid):
    # Define the grid
    grid = [
        ['i', 'i', 'i', 'i', 'i', 'i'],
        ['i', '/', 'i', 'U', 'U', 'i'],
        ['i', 'U', 'U', 'i', 'i', 'i'],
        ['i', 'U', 'U', 'U', 'i', 'i'],
        ['i', 'U', 'U', 'i', 'U', 'i'],
        ['i', 'i', 'i', 'i', 'i', 'i']
    ]
    
    # Directions for moving up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    # Find the start and goal positions
    start = (0, 0)
    goal = None
    for r in range(len(grid)):
        for c in range(len(grid[0])):
            if grid[r][c] == '/':
                goal = (r, c)
                break
        if goal:
            break
    
    # BFS initialization
    queue = deque([(start, 0)])  # (position, steps)
    visited = set([start])
    
    # BFS loop
    while queue:
        (r, c), steps = queue.popleft()
        
        # Check if we reached the goal
        if (r, c) == goal:
            return steps
        
        # Explore neighbors
        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if 0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and (nr, nc) not in visited:
                if grid[nr][nc] == 'U' or grid[nr][nc] == '/':
                    visited.add((nr, nc))
                    queue.append(((nr, nc), steps + 1))
    
    return -1  # If no path is found

# Execute the function and print the result
print(min_steps_to_goal([]))