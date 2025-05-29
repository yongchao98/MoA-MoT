from collections import deque

def find_shortest_path(grid):
    # Directions for moving in the grid: (row_change, col_change, direction)
    directions = [(-1, 0, 'up'), (1, 0, 'down'), (0, -1, 'left'), (0, 1, 'right')]
    
    # Find the start (*) and destination (#) positions
    start = None
    destination = None
    for r in range(len(grid)):
        for c in range(len(grid[0])):
            if grid[r][c] == '*':
                start = (r, c)
            elif grid[r][c] == '#':
                destination = (r, c)
    
    # Initialize BFS
    queue = deque([(start, [])])
    visited = set()
    visited.add(start)
    
    # Perform BFS
    while queue:
        (current_r, current_c), path = queue.popleft()
        
        # Check if we reached the destination
        if (current_r, current_c) == destination:
            return ' '.join(path)
        
        # Explore neighbors
        for dr, dc, direction in directions:
            new_r, new_c = current_r + dr, current_c + dc
            if (0 <= new_r < len(grid) and 0 <= new_c < len(grid[0]) and
                grid[new_r][new_c] != 'X' and (new_r, new_c) not in visited):
                visited.add((new_r, new_c))
                queue.append(((new_r, new_c), path + [direction]))
    
    # If we exhaust the queue without finding the destination
    return "infeasible"

# Define the grid
grid = [
    ['X', 'O', 'O', 'O', '#', 'O', 'O', 'X'],
    ['O', 'O', 'X', 'X', 'O', 'X', 'X', 'O'],
    ['O', 'X', 'O', 'X', 'O', 'O', 'O', 'O'],
    ['O', 'O', 'X', 'O', 'O', 'X', 'O', 'O'],
    ['O', 'O', 'O', 'X', 'X', '*', 'O', 'O']
]

# Find the shortest path
result = find_shortest_path(grid)
print(result)