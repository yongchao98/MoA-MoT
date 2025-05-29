from collections import deque

def find_shortest_path(grid):
    # Directions for moving in the grid: right, down, left, up
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    direction_names = ["right", "down", "left", "up"]
    
    # Find the start (*) and destination (#) positions
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
    
    # BFS setup
    queue = deque([(start, [])])  # (current_position, path_taken)
    visited = set()
    visited.add(start)
    
    while queue:
        (current, path) = queue.popleft()
        
        if current == destination:
            return ' '.join(path)
        
        for direction, (di, dj) in zip(direction_names, directions):
            ni, nj = current[0] + di, current[1] + dj
            if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and grid[ni][nj] in ('O', '#') and (ni, nj) not in visited:
                visited.add((ni, nj))
                queue.append(((ni, nj), path + [direction]))
    
    return "infeasible"

# Define the grid
grid = [
    ['O', 'O', 'O', 'O', 'X', 'O', 'X'],
    ['O', 'O', 'O', 'O', 'X', 'O', 'X'],
    ['X', 'O', 'O', 'X', 'O', 'O', 'X'],
    ['O', 'X', 'O', '#', '*', 'X', 'X'],
    ['O', 'X', 'X', 'X', 'O', 'X', 'O']
]

# Find the shortest path
result = find_shortest_path(grid)
print(result)