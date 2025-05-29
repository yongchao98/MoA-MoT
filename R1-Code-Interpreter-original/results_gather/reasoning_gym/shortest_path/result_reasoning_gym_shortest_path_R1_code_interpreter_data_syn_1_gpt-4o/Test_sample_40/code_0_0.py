from collections import deque

def find_shortest_path(grid):
    # Directions for moving in the grid: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    direction_names = ['up', 'down', 'left', 'right']
    
    # Find the starting and destination points
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
    
    # Initialize the BFS queue
    queue = deque([(start, [])])  # (current_position, path)
    visited = set()
    visited.add(start)
    
    # Perform BFS
    while queue:
        (current, path) = queue.popleft()
        
        # Check if we have reached the destination
        if current == destination:
            return ' '.join(path)
        
        # Explore neighbors
        for direction, (di, dj) in zip(direction_names, directions):
            ni, nj = current[0] + di, current[1] + dj
            if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and grid[ni][nj] != 'X' and (ni, nj) not in visited:
                visited.add((ni, nj))
                queue.append(((ni, nj), path + [direction]))
    
    return "infeasible"

# Define the grid
grid = [
    ['X', 'O', 'X', 'O', 'O', 'O'],
    ['O', 'O', 'X', '#', 'O', 'O'],
    ['X', 'O', 'O', 'O', 'X', 'X'],
    ['X', 'O', 'X', 'X', 'X', 'O'],
    ['O', 'X', '*', 'X', 'O', 'O']
]

# Find the shortest path
result = find_shortest_path(grid)
print(result)