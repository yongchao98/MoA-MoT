from collections import deque

def find_shortest_path(grid):
    # Directions for moving in the grid: right, down, left, up
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    direction_names = ["right", "down", "left", "up"]
    
    # Find the start (*) and destination (#) positions
    start = None
    destination = None
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == '*':
                start = (i, j)
            elif grid[i][j] == '#':
                destination = (i, j)
    
    if not start or not destination:
        return "infeasible"
    
    # BFS initialization
    queue = deque([(start, [])])  # (current_position, path_taken)
    visited = set()
    visited.add(start)
    
    while queue:
        current, path = queue.popleft()
        
        # Check if we reached the destination
        if current == destination:
            return " ".join(path)
        
        # Explore neighbors
        for idx, (di, dj) in enumerate(directions):
            ni, nj = current[0] + di, current[1] + dj
            if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and grid[ni][nj] != 'X' and (ni, nj) not in visited:
                visited.add((ni, nj))
                queue.append(((ni, nj), path + [direction_names[idx]]))
    
    return "infeasible"

# Define the grid
grid = [
    ['X', 'O', 'X', 'X', 'X', 'X'],
    ['X', 'X', 'O', 'O', 'O', 'O'],
    ['O', 'O', 'O', 'O', 'O', 'X'],
    ['X', '#', '*', 'X', 'X', 'O'],
    ['O', 'O', 'O', 'O', 'O', 'X']
]

# Find the shortest path
result = find_shortest_path(grid)
print(result)