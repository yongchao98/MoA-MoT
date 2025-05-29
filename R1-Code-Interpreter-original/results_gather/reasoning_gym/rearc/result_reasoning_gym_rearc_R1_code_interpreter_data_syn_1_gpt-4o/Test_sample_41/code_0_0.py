def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to output grid
    for i, row in enumerate(input_grid):
        # Find the first occurrence of a number that is not the majority in the row
        majority = max(set(row), key=row.count)
        for j, num in enumerate(row):
            if num != majority:
                output_grid[i][j] = 2
                break
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9],
    [4, 9, 4, 9, 9, 9],
    [9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))