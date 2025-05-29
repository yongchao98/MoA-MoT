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
    queue = deque([(start, [])])
    visited = set()
    visited.add(start)
    
    # Perform BFS
    while queue:
        (current, path) = queue.popleft()
        
        # Check if we reached the destination
        if current == destination:
            return ' '.join(path)
        
        # Explore neighbors
        for direction, (di, dj) in directions.items():
            ni, nj = current[0] + di, current[1] + dj
            if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and grid[ni][nj] in ('O', '#') and (ni, nj) not in visited:
                visited.add((ni, nj))
                queue.append(((ni, nj), path + [direction]))
    
    return "infeasible"

# Define the grid
grid = [
    ['O', 'O', 'O', 'X', 'O', 'X', 'O', 'O'],
    ['O', 'X', 'X', '*', 'X', 'X', 'X', 'X'],
    ['X', 'X', 'O', 'O', 'O', 'O', 'O', 'X'],
    ['X', 'O', 'O', 'X', 'O', 'O', 'X', 'X'],
    ['O', 'O', 'O', 'O', 'O', '#', 'X', 'O'],
    ['X', 'O', 'O', 'O', 'O', 'X', 'X', 'O']
]

# Find and print the shortest path
print(find_shortest_path(grid))