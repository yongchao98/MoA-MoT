from collections import deque

def find_shortest_path(grid):
    # Directions for moving in the grid: right, down, left, up
    directions = [(0, 1, 'right'), (1, 0, 'down'), (0, -1, 'left'), (-1, 0, 'up')]
    
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
    visited = set([start])
    
    # Perform BFS
    while queue:
        (current_x, current_y), path = queue.popleft()
        
        # Check if we have reached the destination
        if (current_x, current_y) == destination:
            return ' '.join(path)
        
        # Explore neighbors
        for dx, dy, move in directions:
            new_x, new_y = current_x + dx, current_y + dy
            if (0 <= new_x < len(grid) and 0 <= new_y < len(grid[0]) and
                (new_x, new_y) not in visited and grid[new_x][new_y] in ('O', '#')):
                visited.add((new_x, new_y))
                queue.append(((new_x, new_y), path + [move]))
    
    return "infeasible"

# Define the grid
grid = [
    ['O', 'X', 'O', 'O', 'O', 'O'],
    ['X', 'O', 'O', 'O', 'O', 'O'],
    ['O', 'O', '#', '*', 'O', 'O'],
    ['O', 'X', 'X', 'O', 'X', 'O'],
    ['O', 'O', 'X', 'X', 'O', 'X']
]

# Find the shortest path
result = find_shortest_path(grid)
print(result)