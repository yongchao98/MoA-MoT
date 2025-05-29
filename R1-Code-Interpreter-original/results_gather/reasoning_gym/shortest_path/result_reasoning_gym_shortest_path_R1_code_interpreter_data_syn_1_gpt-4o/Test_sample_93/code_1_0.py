from collections import deque

def find_shortest_path(grid):
    # Find the start (*) and destination (#) positions
    start = None
    destination = None
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == '*':
                start = (i, j)
            elif grid[i][j] == '#':
                destination = (i, j)
    
    # Directions: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    direction_names = ['up', 'down', 'left', 'right']
    
    # BFS setup
    queue = deque([(start, [])])
    visited = set()
    visited.add(start)
    
    # BFS loop
    while queue:
        (current, path) = queue.popleft()
        
        # Check if we reached the destination
        if current == destination:
            return ' '.join(path)
        
        # Explore neighbors
        for d, name in zip(directions, direction_names):
            neighbor = (current[0] + d[0], current[1] + d[1])
            
            # Check if the neighbor is within bounds and not visited
            if (0 <= neighbor[0] < len(grid) and
                0 <= neighbor[1] < len(grid[0]) and
                neighbor not in visited and
                grid[neighbor[0]][neighbor[1]] in ('O', '#')):
                
                visited.add(neighbor)
                queue.append((neighbor, path + [name]))
    
    return "infeasible"

# Define the grid
grid = [
    ['O', 'O', 'O', 'X', 'X', '*'],
    ['O', 'X', 'O', 'O', 'O', 'O'],
    ['X', 'X', 'X', 'O', '#', 'O'],
    ['O', 'X', 'O', 'O', 'O', 'X'],
    ['O', 'O', 'O', 'O', 'X', 'O'],
    ['O', 'O', 'X', 'O', 'X', 'X']
]

# Find the shortest path
result = find_shortest_path(grid)
print(result)