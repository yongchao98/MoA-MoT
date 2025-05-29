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
    queue = deque([(start, [])])  # (position, path)
    visited = set()
    visited.add(start)
    
    # Perform BFS
    while queue:
        (x, y), path = queue.popleft()
        
        # Check if we reached the destination
        if (x, y) == destination:
            return ' '.join(path)
        
        # Explore neighbors
        for direction, (dx, dy) in directions.items():
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] in ('O', '#') and (nx, ny) not in visited:
                visited.add((nx, ny))
                queue.append(((nx, ny), path + [direction]))
    
    return "infeasible"

# Define the grid
grid = [
    ['X', 'O', 'X', 'O', 'X', 'O', 'X'],
    ['O', 'X', 'X', 'X', 'X', 'X', 'O'],
    ['*', 'O', 'X', 'X', 'O', 'O', 'X'],
    ['X', 'O', 'X', 'O', 'O', 'X', '#'],
    ['O', 'X', 'O', 'X', 'O', 'O', 'O']
]

# Find and print the shortest path
print(find_shortest_path(grid))