from collections import deque

def find_shortest_path(grid):
    # Locate the start (*) and destination (#) points
    start = None
    destination = None
    for i, row in enumerate(grid):
        for j, cell in enumerate(row):
            if cell == '*':
                start = (i, j)
            elif cell == '#':
                destination = (i, j)
    
    if not start or not destination:
        return "infeasible"
    
    # Define the possible directions (up, down, left, right)
    directions = {
        'up': (-1, 0),
        'down': (1, 0),
        'left': (0, -1),
        'right': (0, 1)
    }
    
    # Initialize the BFS queue
    queue = deque([(start, [])])
    visited = set()
    visited.add(start)
    
    # Perform BFS
    while queue:
        (current, path) = queue.popleft()
        
        if current == destination:
            return ' '.join(path)
        
        for direction, (di, dj) in directions.items():
            ni, nj = current[0] + di, current[1] + dj
            if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and grid[ni][nj] != 'X' and (ni, nj) not in visited:
                visited.add((ni, nj))
                queue.append(((ni, nj), path + [direction]))
    
    return "infeasible"

# Define the grid
grid = [
    ['O', 'X', 'O', 'O', 'O', 'O', 'X', '*'],
    ['O', 'X', 'X', 'O', 'O', 'O', 'X', 'X'],
    ['X', 'O', 'X', 'O', 'X', 'O', 'X', 'O'],
    ['X', 'X', 'O', 'O', 'O', 'O', 'O', 'O'],
    ['X', 'X', 'O', 'O', 'X', 'X', 'O', 'X'],
    ['X', 'O', 'O', 'X', 'X', 'O', '#', 'O']
]

# Find and print the shortest path
print(find_shortest_path(grid))