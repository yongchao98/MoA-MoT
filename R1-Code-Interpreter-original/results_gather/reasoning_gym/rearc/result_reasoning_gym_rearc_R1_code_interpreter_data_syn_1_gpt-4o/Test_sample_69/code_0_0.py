def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to output grid

    # Define the region to transform based on the pattern observed
    # For simplicity, let's assume the region is defined by the presence of `8`s
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            if input_grid[i][j] == 8:
                # Transform the region around `8`s
                for k in range(max(0, i-1), min(len(input_grid), i+2)):
                    for l in range(max(0, j-1), min(len(input_grid[i]), j+2)):
                        if input_grid[k][l] != 0:  # Keep `0`s unchanged
                            output_grid[k][l] = 8

    return output_grid

# Test input grid
input_grid = [
    [7, 5, 5, 5, 5, 7, 5, 5, 5, 7, 7, 7],
    [7, 0, 7, 5, 5, 7, 7, 5, 7, 7, 5, 5],
    [5, 7, 5, 7, 0, 7, 7, 5, 5, 7, 5, 5],
    [7, 5, 5, 7, 7, 5, 7, 5, 0, 5, 5, 7],
    [5, 8, 8, 8, 8, 5, 0, 5, 5, 5, 0, 7],
    [5, 8, 0, 0, 8, 7, 5, 0, 5, 7, 5, 5],
    [5, 8, 0, 0, 8, 7, 7, 5, 7, 7, 7, 5],
    [7, 8, 0, 0, 8, 5, 5, 7, 7, 7, 5, 7],
    [5, 8, 8, 8, 8, 5, 7, 7, 7, 0, 0, 7],
    [5, 5, 5, 5, 5, 7, 5, 5, 5, 0, 0, 5],
    [5, 5, 5, 0, 5, 7, 7, 7, 7, 0, 0, 5],
    [5, 5, 7, 5, 5, 5, 5, 7, 5, 7, 5, 7]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))