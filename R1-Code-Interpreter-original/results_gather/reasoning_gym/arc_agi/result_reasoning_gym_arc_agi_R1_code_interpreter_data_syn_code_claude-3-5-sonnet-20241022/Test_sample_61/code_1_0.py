import numpy as np
from collections import deque

def get_connected_twos(grid, start_x, start_y):
    rows, cols = len(grid), len(grid[0])
    region = set()
    queue = deque([(start_x, start_y)])
    
    while queue:
        x, y = queue.popleft()
        if (x, y) in region:
            continue
        if x < 0 or x >= rows or y < 0 or y >= cols:
            continue
        if grid[x][y] != 2:
            continue
            
        region.add((x, y))
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            queue.append((x + dx, y + dy))
    
    return region

def expand_region(grid, output, start_points, value):
    rows, cols = len(grid), len(grid[0])
    queue = deque(list(start_points))
    visited = set(start_points)
    
    while queue:
        x, y = queue.popleft()
        output[x][y] = value
        
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = x + dx, y + dy
            if ((nx, ny) not in visited and 
                0 <= nx < rows and 0 <= ny < cols and 
                grid[nx][ny] in [2, 7] and 
                output[nx][ny] == 9):
                visited.add((nx, ny))
                queue.append((nx, ny))

def solve_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = np.full((rows, cols), 9)
    
    # First find all special numbers
    special_numbers = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] not in [2, 7]:
                special_numbers.append((grid[i][j], i, j))
    
    # Sort special numbers by value
    special_numbers.sort()
    
    # Process each special number and its region
    for value, x, y in special_numbers:
        expand_region(grid, output, {(x, y)}, value)
    
    # Find and process 2's regions
    visited_twos = set()
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 2 and (i, j) not in visited_twos:
                two_region = get_connected_twos(grid, i, j)
                visited_twos.update(two_region)
                
                # Find nearest special number
                min_dist = float('inf')
                nearest_value = 0
                
                for x, y in two_region:
                    for val, sx, sy in special_numbers:
                        dist = abs(x - sx) + abs(y - sy)
                        if dist < min_dist:
                            min_dist = dist
                            nearest_value = val
                
                # Expand the region with the nearest value
                expand_region(grid, output, two_region, nearest_value)
    
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