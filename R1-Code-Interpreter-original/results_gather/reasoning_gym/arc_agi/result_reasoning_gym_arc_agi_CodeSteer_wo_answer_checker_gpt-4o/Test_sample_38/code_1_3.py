def shift_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid
    for i, row in enumerate(input_grid):
        if 9 in row:
            # Find the first occurrence of 9
            first_nine_index = row.index(9)
            # Determine the shift amount based on the first occurrence of 9
            shift_amount = first_nine_index
            # Shift the row to the right by the shift amount
            new_row = [0] * len(row)
            for j, value in enumerate(row):
                if value != 0:
                    new_position = j + shift_amount
                    if new_position < len(row):
                        new_row[new_position] = value
            output_grid[i] = new_row
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 9, 9, 9, 0],
    [0, 0, 0, 9, 6, 9, 0],
    [0, 0, 0, 9, 9, 9, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0]
]

output_grid = shift_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))