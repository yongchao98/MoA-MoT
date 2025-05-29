def transform_grid(input_grid):
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]

    # Get the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])

    # Iterate over each cell in the grid
    for i in range(rows):
        for j in range(cols):
            # Check if the current cell is a non-2 number and not 5
            if input_grid[i][j] != 2 and input_grid[i][j] != 5:
                # Apply the cross pattern
                num = input_grid[i][j]
                output_grid[i][j] = num
                if i > 0:
                    output_grid[i-1][j] = num
                if i < rows - 1:
                    output_grid[i+1][j] = num
                if j > 0:
                    output_grid[i][j-1] = num
                if j < cols - 1:
                    output_grid[i][j+1] = num

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

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))