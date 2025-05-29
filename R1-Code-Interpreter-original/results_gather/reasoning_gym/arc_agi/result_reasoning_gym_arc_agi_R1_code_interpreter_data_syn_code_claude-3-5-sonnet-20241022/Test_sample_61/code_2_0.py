import numpy as np
from collections import deque

def flood_fill(grid, x, y, old_value, new_value, output):
    if old_value == new_value:
        return
    
    rows, cols = len(grid), len(grid[0])
    queue = deque([(x, y)])
    
    while queue:
        cx, cy = queue.popleft()
        if (cx < 0 or cx >= rows or cy < 0 or cy >= cols or 
            grid[cx][cy] != old_value or output[cx][cy] == new_value):
            continue
            
        output[cx][cy] = new_value
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = cx + dx, cy + dy
            queue.append((nx, ny))

def solve_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = np.full((rows, cols), -1)
    
    # First pass: identify special numbers and their regions
    special_positions = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 7 and grid[i][j] != 2:
                special_positions.append((i, j, grid[i][j]))
    
    # Sort special positions by value
    special_positions.sort(key=lambda x: x[2])
    
    # Process each special number
    for x, y, value in special_positions:
        if output[x][y] == -1:
            flood_fill(grid, x, y, 7, value, output)
            
    # Process regions marked by 2's
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 2:
                # Find adjacent special numbers
                adjacent_value = None
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = i + dx, j + dy
                    if (0 <= nx < rows and 0 <= ny < cols and 
                        grid[nx][ny] not in [2, 7] and output[nx][ny] != -1):
                        adjacent_value = output[nx][ny]
                        break
                
                if adjacent_value is None:
                    adjacent_value = 0  # Default value if no special numbers found
                
                flood_fill(grid, i, j, 2, adjacent_value, output)
    
    # Fill remaining areas (where value is still -1) with 9
    for i in range(rows):
        for j in range(cols):
            if output[i][j] == -1:
                output[i][j] = 9
    
    return '\n'.join(' '.join(map(str, row)) for row in output)

# Test input
test_input = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 7, 6, 7],
    [7, 7, 2, 2, 2, 2, 2, 7, 7, 2, 7, 7, 7],
    [7, 7, 2, 7, 7, 7, 2, 7, 7, 2, 2, 7, 7],
    [7, 2, 2, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2],
    [7, 2, 7, 7, 7, 2, 7, 7, 7, 7, 5, 7, 7],
    [7, 2, 7, 9, 7, 2, 7, 7, 7, 7, 7, 7, 7],
    [7, 2, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7],
    [7, 8, 7, 7, 7, 2, 2, 2, 7, 7, 7, 7, 7],
    [7, 2, 7, 7, 7, 7, 7, 2, 7, 7, 7, 7, 7],
    [7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [7, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 2, 2, 7, 7, 7, 3, 7, 7, 7],
    [7, 7, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7]
]

print(solve_grid(test_input))