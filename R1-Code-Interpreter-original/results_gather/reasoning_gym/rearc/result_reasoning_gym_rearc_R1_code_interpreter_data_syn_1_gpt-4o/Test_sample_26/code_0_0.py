def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    for i in range(1, len(input_grid)):
        for j in range(len(input_grid[i])):
            # Check if the current row has a sequence of 4s
            if input_grid[i][j] == 4:
                # Replace the sequence with the pattern from the row above
                output_grid[i] = input_grid[i - 1]
                break  # No need to check further in this row

    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 1, 6, 7, 1, 6, 7, 1, 6, 7],
    [7, 6, 1, 7, 6, 1, 7, 6, 1, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 1, 6, 7, 1, 6, 7, 1, 6, 7],
    [7, 6, 1, 7, 6, 1, 7, 6, 1, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 1, 6, 7, 1, 6, 7, 1, 6, 7],
    [7, 6, 4, 4, 4, 1, 7, 6, 1, 7],
    [7, 7, 4, 4, 4, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))