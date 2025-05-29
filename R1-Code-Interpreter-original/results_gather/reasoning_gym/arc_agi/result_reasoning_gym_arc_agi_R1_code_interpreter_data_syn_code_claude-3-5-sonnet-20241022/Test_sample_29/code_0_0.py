import numpy as np
from collections import deque

def find_connected_components(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    components = []
    
    def bfs(r, c, target):
        component = set()
        queue = deque([(r, c)])
        component.add((r, c))
        
        while queue:
            r, c = queue.popleft()
            for nr, nc in [(r+1,c), (r-1,c), (r,c+1), (r,c-1)]:
                if (0 <= nr < rows and 0 <= nc < cols and 
                    grid[nr][nc] == target and 
                    (nr, nc) not in component):
                    queue.append((nr, nc))
                    component.add((nr, nc))
        return component

    for i in range(rows):
        for j in range(cols):
            if (i, j) not in visited and grid[i][j] != 8:
                component = bfs(i, j, grid[i][j])
                visited.update(component)
                components.append((grid[i][j], component))
    
    return components

def fill_output_grid(input_grid, components):
    rows, cols = len(input_grid), len(input_grid[0])
    output = [[3 for _ in range(cols)] for _ in range(rows)]  # Default background
    
    # Fill components
    for value, positions in components:
        fill_value = value
        for r, c in positions:
            output[r][c] = fill_value
            
        # Expand region
        queue = deque(positions)
        visited = set(positions)
        while queue:
            r, c = queue.popleft()
            for nr, nc in [(r+1,c), (r-1,c), (r,c+1), (r,c-1)]:
                if (0 <= nr < rows and 0 <= nc < cols and 
                    (nr, nc) not in visited and 
                    input_grid[nr][nc] == 8):
                    output[nr][nc] = fill_value
                    queue.append((nr, nc))
                    visited.add((nr, nc))
    
    return output

# Input grid
input_grid = [
    [8,8,8,8,8,8,8,2,8,8,8,8,8],
    [8,8,8,0,8,8,8,2,2,8,8,8,8],
    [8,8,8,8,8,8,8,8,2,8,8,8,8],
    [2,2,2,2,2,2,2,2,2,2,2,2,8],
    [8,8,8,8,8,2,8,8,8,8,8,2,8],
    [8,8,8,8,8,2,2,2,8,8,8,6,8],
    [8,8,8,8,8,8,8,2,8,8,8,2,8],
    [8,8,8,8,8,8,8,2,8,4,8,2,8],
    [8,8,9,8,8,8,8,2,8,8,8,2,8],
    [2,2,2,8,8,8,2,2,8,8,2,2,8],
    [8,8,2,2,8,8,2,8,8,8,2,8,8],
    [8,8,8,2,8,8,2,2,2,2,2,8,8],
    [8,7,8,3,8,8,8,8,8,8,8,8,8]
]

# Process the grid
components = find_connected_components(input_grid)
output_grid = fill_output_grid(input_grid, components)

# Print the result
for row in output_grid:
    print(' '.join(map(str, row)))