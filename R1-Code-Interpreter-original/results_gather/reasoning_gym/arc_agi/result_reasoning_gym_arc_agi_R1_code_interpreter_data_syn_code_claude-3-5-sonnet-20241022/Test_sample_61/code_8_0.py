import numpy as np
from collections import deque

def get_connected_region(grid, start_x, start_y):
    rows, cols = len(grid), len(grid[0])
    value = grid[start_x][start_y]
    region = set()
    queue = deque([(start_x, start_y)])
    
    while queue:
        x, y = queue.popleft()
        if (x, y) in region:
            continue
        if x < 0 or x >= rows or y < 0 or y >= cols:
            continue
        if grid[x][y] != value:
            continue
            
        region.add((x, y))
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            queue.append((x + dx, y + dy))
    
    return region

def flood_fill(output, start_x, start_y, value, boundary):
    rows, cols = len(output), len(output[0])
    if (start_x, start_y) in boundary:
        return
        
    queue = deque([(start_x, start_y)])
    filled = set()
    
    while queue:
        x, y = queue.popleft()
        if (x, y) in filled:
            continue
        if x < 0 or x >= rows or y < 0 or y >= cols:
            continue
        if (x, y) in boundary:
            continue
            
        output[x][y] = value
        filled.add((x, y))
        
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = x + dx, y + dy
            if (nx, ny) not in boundary:
                queue.append((nx, ny))

def solve_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = np.full((rows, cols), 9)
    boundary = set()
    
    # First handle special numbers (non-2 and non-7)
    special_regions = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] not in [2, 7]:
                special_regions.append((i, j, grid[i][j]))
    
    # Process 2's regions
    two_regions = []
    visited_twos = set()
    
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 2 and (i, j) not in visited_twos:
                region = get_connected_region(grid, i, j)
                visited_twos.update(region)
                
                # Find nearest special number
                min_dist = float('inf')
                nearest_value = 0
                
                for x, y in region:
                    for sx, sy, val in special_regions:
                        dist = abs(x - sx) + abs(y - sy)
                        if dist < min_dist:
                            min_dist = dist
                            nearest_value = val
                
                two_regions.append((region, nearest_value))
                boundary.update(region)
    
    # Fill special numbers first
    for x, y, val in special_regions:
        output[x][y] = val
        boundary.add((x, y))
        flood_fill(output, x, y, val, boundary)
    
    # Fill 2's regions
    for region, value in two_regions:
        for x, y in region:
            output[x][y] = value
        for x, y in region:
            flood_fill(output, x, y, value, boundary)
    
    # Fill remaining areas with 9
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in boundary and output[i][j] != 9:
                flood_fill(output, i, j, 9, boundary)
    
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