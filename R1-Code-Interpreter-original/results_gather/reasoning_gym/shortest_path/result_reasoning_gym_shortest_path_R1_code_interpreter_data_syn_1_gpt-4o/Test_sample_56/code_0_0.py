from collections import deque

def find_shortest_path(grid):
    # Directions for moving in the grid
    directions = {'up': (-1, 0), 'down': (1, 0), 'left': (0, -1), 'right': (0, 1)}
    
    # Find the start (*) and destination (#) positions
    start = None
    destination = None
    for i, row in enumerate(grid):
        for j, cell in enumerate(row):
            if cell == '*':
                start = (i, j)
            elif cell == '#':
                destination = (i, j)
    
    # Initialize BFS
    queue = deque([(start, [])])  # (current_position, path_taken)
    visited = set()
    visited.add(start)
    
    # Perform BFS
    while queue:
        (current_x, current_y), path = queue.popleft()
        
        # Check if we reached the destination
        if (current_x, current_y) == destination:
            return ' '.join(path)
        
        # Explore neighbors
        for direction, (dx, dy) in directions.items():
            new_x, new_y = current_x + dx, current_y + dy
            if 0 <= new_x < len(grid) and 0 <= new_y < len(grid[0]) and grid[new_x][new_y] in ('O', '#') and (new_x, new_y) not in visited:
                visited.add((new_x, new_y))
                queue.append(((new_x, new_y), path + [direction]))
    
    # If we exhaust the queue without finding the destination
    return "infeasible"

# Define the grid
grid = [
    ['X', 'X', 'X', 'O', 'X', 'X'],
    ['O', 'X', 'O', 'O', 'O', 'O'],
    ['X', 'O', 'X', 'X', 'O', 'X'],
    ['X', 'O', 'O', 'X', 'X', 'X'],
    ['X', 'O', 'O', 'O', 'O', 'O'],
    ['O', 'O', 'O', 'X', '#', 'X'],
    ['X', 'O', 'O', 'O', 'O', 'X'],
    ['O', 'X', 'X', 'X', 'O', '*']
]

# Find and print the shortest path
print(find_shortest_path(grid))