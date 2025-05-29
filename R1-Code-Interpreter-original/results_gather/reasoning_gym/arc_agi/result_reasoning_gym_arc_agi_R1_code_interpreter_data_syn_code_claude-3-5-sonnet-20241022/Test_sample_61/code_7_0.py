import numpy as np
from collections import deque

def get_regions(grid, value):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    regions = []
    
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == value and (i, j) not in visited:
                region = set()
                queue = deque([(i, j)])
                
                while queue:
                    x, y = queue.popleft()
                    if (x, y) in visited:
                        continue
                    if x < 0 or x >= rows or y < 0 or y >= cols:
                        continue
                    if grid[x][y] != value:
                        continue
                    
                    visited.add((x, y))
                    region.add((x, y))
                    
                    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                        queue.append((x + dx, y + dy))
                
                if region:
                    regions.append(region)
    return regions

def find_nearest_special(grid, pos, special_positions):
    x, y = pos
    min_dist = float('inf')
    nearest_value = 0
    
    for sx, sy, val in special_positions:
        dist = abs(x - sx) + abs(y - sy)
        if dist < min_dist:
            min_dist = dist
            nearest_value = val
    
    return nearest_value

def solve_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = np.full((rows, cols), 9)
    
    # Find special numbers
    special_positions = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] not in [2, 7]:
                special_positions.append((i, j, grid[i][j]))
                output[i][j] = grid[i][j]
    
    # Process 2's regions
    two_regions = get_regions(grid, 2)
    for region in two_regions:
        # Find value for this region
        value = find_nearest_special(grid, next(iter(region)), special_positions)
        
        # Fill region and expand
        queue = deque(list(region))
        visited = set(region)
        
        while queue:
            x, y = queue.popleft()
            output[x][y] = value
            
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nx, ny = x + dx, y + dy
                if (0 <= nx < rows and 0 <= ny < cols and 
                    (nx, ny) not in visited and 
                    (grid[nx][ny] == 7 or grid[nx][ny] == 2)):
                    visited.add((nx, ny))
                    queue.append((nx, ny))
                    output[nx][ny] = value
    
    # Expand special number regions
    for sx, sy, val in special_positions:
        queue = deque([(sx, sy)])
        visited = {(sx, sy)}
        
        while queue:
            x, y = queue.popleft()
            output[x][y] = val
            
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nx, ny = x + dx, y + dy
                if (0 <= nx < rows and 0 <= ny < cols and 
                    (nx, ny) not in visited and 
                    output[nx][ny] == 9 and 
                    grid[nx][ny] == 7):
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