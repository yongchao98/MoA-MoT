import numpy as np
from scipy.stats import mode

def apply_rule(input_grid):
    output_grid = np.array(input_grid)
    rows, cols = output_grid.shape
    
    for i in range(rows):
        for j in range(cols):
            # Get the surrounding elements
            surrounding = []
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if di == 0 and dj == 0:
                        continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < rows and 0 <= nj < cols:
                        surrounding.append(output_grid[ni, nj])
            
            # If the current element is different from the mode of its surroundings, replace it
            if surrounding:
                most_common = mode(surrounding).mode[0]
                if output_grid[i, j] != most_common:
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

output_grid = apply_rule(input_grid)
print(output_grid)