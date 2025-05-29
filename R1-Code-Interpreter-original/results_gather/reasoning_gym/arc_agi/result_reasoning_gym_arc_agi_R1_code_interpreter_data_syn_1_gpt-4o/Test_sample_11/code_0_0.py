def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to output grid
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    for i in range(num_rows):
        for j in range(num_cols):
            if input_grid[i][j] == 0:
                # Replace 0 with the number at the start of the row
                output_grid[i][j] = input_grid[i][0]
    
    return output_grid

# Test input grid
input_grid = [
    [4, 8, 8, 8, 8, 8, 8, 8, 8, 6],
    [4, 8, 0, 8, 8, 8, 0, 8, 8, 6],
    [4, 8, 8, 8, 0, 8, 8, 8, 8, 6],
    [4, 8, 8, 8, 8, 8, 8, 8, 8, 6],
    [4, 8, 8, 0, 8, 8, 8, 8, 8, 6],
    [4, 8, 8, 8, 8, 8, 0, 8, 8, 6],
    [4, 8, 8, 0, 8, 8, 8, 8, 8, 6],
    [4, 8, 8, 8, 8, 8, 8, 8, 8, 6],
    [4, 8, 8, 8, 8, 0, 8, 8, 0, 6],
    [4, 0, 8, 8, 8, 8, 8, 8, 8, 6]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))