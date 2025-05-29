from collections import deque

def find_shortest_path(grid):
    # Define the grid size
    rows, cols = len(grid), len(grid[0])
    
    # Locate the start (*) and destination (#)
    start = None
    destination = None
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '*':
                start = (r, c)
            elif grid[r][c] == '#':
                destination = (r, c)
    
    # Directions: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    direction_names = ['up', 'down', 'left', 'right']
    
    # BFS setup
    queue = deque([(start, [])])  # (current_position, path_taken)
    visited = set()
    visited.add(start)
    
    while queue:
        (current_r, current_c), path = queue.popleft()
        
        # Check if we reached the destination
        if (current_r, current_c) == destination:
            return ' '.join(path)
        
        # Explore neighbors
        for i, (dr, dc) in enumerate(directions):
            new_r, new_c = current_r + dr, current_c + dc
            if 0 <= new_r < rows and 0 <= new_c < cols and grid[new_r][new_c] != 'X' and (new_r, new_c) not in visited:
                visited.add((new_r, new_c))
                queue.append(((new_r, new_c), path + [direction_names[i]]))
    
    return "infeasible"

# Define the grid
grid = [
    ['O', 'O', 'O', 'X', 'X', '*'],
    ['O', 'X', 'O', 'O', 'O', 'O'],
    ['X', 'X', 'X', 'O', '#', 'O'],
    ['O', 'X', 'O', 'O', 'O', 'X'],
    ['O', 'O', 'O', 'O', 'X', 'O'],
    ['O', 'O', 'X', 'O', 'X', 'X']
]

# Find and print the shortest path
print(find_shortest_path(grid))