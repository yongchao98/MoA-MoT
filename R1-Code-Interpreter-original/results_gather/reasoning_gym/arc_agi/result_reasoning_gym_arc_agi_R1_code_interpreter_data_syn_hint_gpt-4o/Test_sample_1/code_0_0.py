def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        # Create two rows for each input row
        new_row1 = row + row[::-1]  # Original row + reversed row
        new_row2 = row[::-1] + row  # Reversed row + original row
        output_grid.append(new_row1)
        output_grid.append(new_row2)
    return output_grid

# Test input grid
input_grid = [
    [4, 6, 4, 2],
    [2, 5, 2, 4],
    [6, 5, 6, 6],
    [6, 5, 6, 6]
]

# Generate the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))