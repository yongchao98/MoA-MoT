import numpy as np
from collections import deque

def get_distance(p1, p2):
    return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])

def solve_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = np.full((rows, cols), -1)
    
    # First, find all special numbers and their positions
    special_positions = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] not in [2, 7]:
                special_positions.append((i, j, grid[i][j]))
    
    # Find all connected regions of 2's
    two_regions = []
    visited = set()
    
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 2 and (i, j) not in visited:
                region = set()
                queue = deque([(i, j)])
                
                while queue:
                    x, y = queue.popleft()
                    if (x, y) in visited:
                        continue
                    if x < 0 or x >= rows or y < 0 or y >= cols:
                        continue
                    if grid[x][y] != 2:
                        continue
                        
                    visited.add((x, y))
                    region.add((x, y))
                    
                    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                        queue.append((x + dx, y + dy))
                
                if region:
                    two_regions.append(region)
    
    # Process each 2's region
    for region in two_regions:
        # Find nearest special number
        min_dist = float('inf')
        fill_value = 0  # Default value
        
        for x, y in region:
            for sx, sy, val in special_positions:
                dist = get_distance((x, y), (sx, sy))
                if dist < min_dist:
                    min_dist = dist
                    fill_value = val
        
        # Fill the region
        for x, y in region:
            output[x][y] = fill_value
    
    # Fill special numbers
    for x, y, val in special_positions:
        output[x][y] = val
    
    # Fill remaining spaces (background)
    for i in range(rows):
        for j in range(cols):
            if output[i][j] == -1:
                # Find nearest filled position
                min_dist = float('inf')
                fill_value = 9  # Default background
                
                for x in range(rows):
                    for y in range(cols):
                        if output[x][y] != -1:
                            dist = get_distance((i, j), (x, y))
                            if dist < min_dist:
                                min_dist = dist
                                fill_value = output[x][y]
                
                output[i][j] = fill_value
    
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