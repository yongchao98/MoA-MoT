import numpy as np
from collections import deque

def get_connected_region(grid, start_x, start_y, target):
    rows, cols = len(grid), len(grid[0])
    region = set()
    queue = deque([(start_x, start_y)])
    
    while queue:
        x, y = queue.popleft()
        if (x, y) in region:
            continue
        if x < 0 or x >= rows or y < 0 or y >= cols:
            continue
        if grid[x][y] != target:
            continue
            
        region.add((x, y))
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            queue.append((x + dx, y + dy))
    
    return region

def find_nearest_special(grid, region):
    special_values = []
    for x, y in region:
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = x + dx, y + dy
            if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and 
                grid[nx][ny] not in [2, 7]):
                special_values.append(grid[nx][ny])
    return special_values[0] if special_values else 0

def solve_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = np.full((rows, cols), 9)  # Initialize with background value
    processed = set()
    
    # First find all special numbers and their regions
    special_positions = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] not in [2, 7]:
                special_positions.append((i, j))
                
    # Process each special number
    for sx, sy in special_positions:
        value = grid[sx][sy]
        output[sx][sy] = value
        
        # Expand the region
        queue = deque([(sx, sy)])
        visited = {(sx, sy)}
        
        while queue:
            x, y = queue.popleft()
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nx, ny = x + dx, y + dy
                if ((nx, ny) not in visited and 0 <= nx < rows and 0 <= ny < cols and 
                    grid[nx][ny] == 7):
                    output[nx][ny] = value
                    visited.add((nx, ny))
                    queue.append((nx, ny))
    
    # Process regions of 2's
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 2 and (i, j) not in processed:
                region = get_connected_region(grid, i, j, 2)
                processed.update(region)
                
                # Find the value to fill this region with
                fill_value = find_nearest_special(grid, region)
                
                # Fill the region and expand
                for rx, ry in region:
                    output[rx][ry] = fill_value
                    
                # Expand into adjacent 7's
                for rx, ry in region:
                    queue = deque([(rx, ry)])
                    visited = {(rx, ry)}
                    
                    while queue:
                        x, y = queue.popleft()
                        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                            nx, ny = x + dx, y + dy
                            if ((nx, ny) not in visited and 0 <= nx < rows and 
                                0 <= ny < cols and grid[nx][ny] == 7 and 
                                output[nx][ny] == 9):
                                output[nx][ny] = fill_value
                                visited.add((nx, ny))
                                queue.append((nx, ny))
    
    return output

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

result = solve_grid(test_input)
print('\n'.join(' '.join(map(str, row)) for row in result))