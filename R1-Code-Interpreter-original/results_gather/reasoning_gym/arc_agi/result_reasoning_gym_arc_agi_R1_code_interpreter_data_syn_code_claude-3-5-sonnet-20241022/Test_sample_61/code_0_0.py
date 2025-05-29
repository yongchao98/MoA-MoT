import numpy as np
from collections import deque

def flood_fill(grid, x, y, target_value):
    if x < 0 or x >= len(grid) or y < 0 or y >= len(grid[0]):
        return set()
    
    filled = set()
    queue = deque([(x, y)])
    while queue:
        curr_x, curr_y = queue.popleft()
        if ((curr_x, curr_y) in filled or 
            curr_x < 0 or curr_x >= len(grid) or 
            curr_y < 0 or curr_y >= len(grid[0]) or
            grid[curr_x][curr_y] != target_value):
            continue
            
        filled.add((curr_x, curr_y))
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            queue.append((curr_x + dx, curr_y + dy))
    return filled

def find_regions(grid):
    regions = []
    visited = set()
    
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if (i, j) not in visited and grid[i][j] != 7:
                region = flood_fill(grid, i, j, grid[i][j])
                if region:
                    regions.append((grid[i][j], region))
                visited.update(region)
    return regions

def solve_grid(input_grid):
    # Convert input to numpy array
    grid = np.array([list(map(int, row)) for row in input_grid.strip().split('\n')])
    output = np.zeros_like(grid)
    
    # Find all regions of non-7 numbers
    regions = find_regions(grid)
    
    # Sort regions by size (larger regions first)
    regions.sort(key=lambda x: len(x[1]), reverse=True)
    
    # Process each region
    for value, region in regions:
        region_number = value
        if value == 2:  # Special case for 2's
            region_number = 0  # Default fill for 2's regions
            # Check nearby special numbers
            for x, y in region:
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and 
                        grid[nx][ny] not in [2, 7]):
                        region_number = grid[nx][ny]
                        break
                if region_number != 0:
                    break
        
        # Fill the region
        for x, y in region:
            output[x][y] = region_number
            
        # Expand region
        expanded = set()
        for x, y in region:
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nx, ny = x + dx, y + dy
                if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and 
                    grid[nx][ny] == 7 and output[nx][ny] == 0):
                    output[nx][ny] = region_number
                    expanded.add((nx, ny))
    
    # Fill remaining areas
    default_fill = 9
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if output[i][j] == 0:
                output[i][j] = default_fill
    
    return '\n'.join(' '.join(map(str, row)) for row in output.astype(int))

# Test input
test_input = """7 7 7 7 7 7 7 7 7 0 7 6 7
7 7 2 2 2 2 2 7 7 2 7 7 7
7 7 2 7 7 7 2 7 7 2 2 7 7
7 2 2 7 7 2 2 7 7 7 2 2 2
7 2 7 7 7 2 7 7 7 7 5 7 7
7 2 7 9 7 2 7 7 7 7 7 7 7
7 2 7 7 7 2 7 7 7 7 7 7 7
7 8 7 7 7 2 2 2 7 7 7 7 7
7 2 7 7 7 7 7 2 7 7 7 7 7
7 2 2 2 2 2 2 2 2 2 2 2 2
7 7 7 7 2 7 7 7 7 7 7 7 7
7 7 7 7 2 2 7 7 7 3 7 7 7
7 7 7 7 7 2 7 7 7 7 7 7 7"""

print(solve_grid(test_input))