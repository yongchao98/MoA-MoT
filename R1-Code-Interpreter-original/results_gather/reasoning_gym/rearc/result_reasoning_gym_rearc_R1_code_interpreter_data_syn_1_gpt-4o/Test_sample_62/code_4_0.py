import numpy as np
from collections import Counter

def replace_outliers(grid):
    grid = np.array(grid)
    rows, cols = grid.shape
    output_grid = grid.copy()
    
    for i in range(1, rows-1):
        for j in range(1, cols-1):
            # Get the surrounding numbers
            surrounding = [
                grid[i-1, j-1], grid[i-1, j], grid[i-1, j+1],
                grid[i, j-1],               grid[i, j+1],
                grid[i+1, j-1], grid[i+1, j], grid[i+1, j+1]
            ]
            # Find the most common surrounding number
            most_common = Counter(surrounding).most_common(1)[0][0]
            # Replace if the current number is different from the most common
            if grid[i, j] != most_common:
                output_grid[i, j] = most_common
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 7, 7],
    [8, 8, 8, 8, 8, 7, 7],
    [8, 8, 8, 8, 8, 7, 1],
    [8, 8, 8, 8, 8, 7, 7],
    [8, 8, 2, 8, 8, 7, 7],
    [8, 8, 5, 8, 8, 7, 7],
    [8, 8, 8, 8, 8, 7, 7]
]

output_grid = replace_outliers(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))