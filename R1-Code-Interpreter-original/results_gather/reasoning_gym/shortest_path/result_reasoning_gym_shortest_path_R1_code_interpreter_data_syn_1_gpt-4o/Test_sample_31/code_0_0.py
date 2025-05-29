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
    
    # BFS initialization
    queue = deque([(start, [])])  # (current_position, path_taken)
    visited = set()
    visited.add(start)
    
    # Perform BFS
    while queue:
        (current, path) = queue.popleft()
        
        # Check if we reached the destination
        if current == destination:
            return "<<<" + " ".join(path) + ">>>"
        
        # Explore neighbors
        for direction, (di, dj) in zip(direction_names, directions):
            ni, nj = current[0] + di, current[1] + dj
            if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and grid[ni][nj] != 'X' and (ni, nj) not in visited:
                visited.add((ni, nj))
                queue.append(((ni, nj), path + [direction]))
    
    return "infeasible"

# Define the grid
grid = [
    ['X', 'O', 'O', 'X', 'O', 'X', 'X', 'X'],
    ['O', 'O', 'O', 'X', 'O', 'O', 'X', 'X'],
    ['O', 'O', '*', 'O', 'O', 'X', 'O', 'O'],
    ['O', 'O', 'X', 'O', 'X', 'X', 'X', 'O'],
    ['#', 'O', 'O', 'O', 'O', 'X', 'O', 'X'],
    ['X', 'X', 'X', 'O', 'O', 'O', 'O', 'X'],
    ['X', 'O', 'O', 'O', 'O', 'O', 'X', 'X'],
    ['O', 'X', 'O', 'X', 'X', 'O', 'O', 'O']
]

# Find and print the shortest path
print(find_shortest_path(grid))