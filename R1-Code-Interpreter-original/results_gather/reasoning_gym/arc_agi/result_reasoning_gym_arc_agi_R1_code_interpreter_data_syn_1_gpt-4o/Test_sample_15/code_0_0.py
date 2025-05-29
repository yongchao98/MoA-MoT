def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 1:
                # Add 1s in a cross pattern around the 1
                if i > 0:
                    output_grid[i-1][j] = 1
                if i < rows - 1:
                    output_grid[i+1][j] = 1
                if j > 0:
                    output_grid[i][j-1] = 1
                if j < cols - 1:
                    output_grid[i][j+1] = 1
    
    return output_grid

# Test input grid
input_grid = [
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 0, 1, 0, 0, 0, 2],
    [2, 2, 2, 2, 1, 0, 1, 2, 2, 2],
    [2, 2, 2, 2, 0, 1, 0, 0, 0, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))