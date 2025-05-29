import numpy as np
from collections import deque

def flood_fill(grid, x, y, value, output, visited):
    rows, cols = len(grid), len(grid[0])
    queue = deque([(x, y)])
    region = set()
    
    while queue:
        cx, cy = queue.popleft()
        if ((cx, cy) in visited or cx < 0 or cx >= rows or cy < 0 or cy >= cols or 
            (grid[cx][cy] != 7 and grid[cx][cy] != 2)):
            continue
            
        visited.add((cx, cy))
        region.add((cx, cy))
        output[cx][cy] = value
        
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = cx + dx, cy + dy
            queue.append((nx, ny))
    return region

def solve_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = np.full((rows, cols), 9)  # Default background
    visited = set()
    
    # First, process special numbers (non-2 and non-7)
    special_numbers = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] not in [2, 7]:
                special_numbers.append((i, j, grid[i][j]))
    
    # Process special numbers
    for x, y, value in special_numbers:
        if (x, y) not in visited:
            output[x][y] = value
            flood_fill(grid, x, y, value, output, visited)
    
    # Process 2's regions
    twos_regions = []
    visited_twos = set()
    
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 2 and (i, j) not in visited_twos:
                region = set()
                queue = deque([(i, j)])
                
                while queue:
                    cx, cy = queue.popleft()
                    if ((cx, cy) in visited_twos or cx < 0 or cx >= rows or 
                        cy < 0 or cy >= cols or grid[cx][cy] != 2):
                        continue
                    
                    visited_twos.add((cx, cy))
                    region.add((cx, cy))
                    
                    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                        nx, ny = cx + dx, cy + dy
                        queue.append((nx, ny))
                
                if region:
                    # Find adjacent special numbers
                    fill_value = 0  # Default fill value for 2's regions
                    for rx, ry in region:
                        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                            nx, ny = rx + dx, ry + dy
                            if (0 <= nx < rows and 0 <= ny < cols and 
                                grid[nx][ny] not in [2, 7]):
                                fill_value = grid[nx][ny]
                                break
                        if fill_value != 0:
                            break
                    
                    # Fill the region and expand
                    for rx, ry in region:
                        output[rx][ry] = fill_value
                    for rx, ry in region:
                        flood_fill(grid, rx, ry, fill_value, output, visited)
    
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