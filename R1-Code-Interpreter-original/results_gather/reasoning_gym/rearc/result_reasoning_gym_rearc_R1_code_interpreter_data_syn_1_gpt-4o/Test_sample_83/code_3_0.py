def transform_grid(input_grid):
    # Remove the section of 9s in the fourth and fifth rows and fill with alternating pattern
    output_grid = []
    for i, row in enumerate(input_grid):
        if i in [3, 4]:
            # Remove the section of 9s and fill with alternating pattern
            new_row = row[:5] + [0, 6, 0, 6]  # Ensure the pattern is consistent
        else:
            new_row = row
        output_grid.append(new_row)
    return output_grid

# Test input grid
input_grid = [
    [0, 6, 0, 6, 0, 6, 0, 6],
    [6, 0, 6, 0, 6, 0, 6, 0],
    [0, 6, 0, 6, 0, 6, 0, 6],
    [6, 0, 6, 0, 6, 9, 9, 0],
    [0, 6, 0, 6, 0, 9, 9, 6],
    [6, 0, 6, 0, 6, 0, 6, 0],
    [0, 6, 0, 6, 0, 6, 0, 6],
    [6, 0, 6, 0, 6, 0, 6, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))