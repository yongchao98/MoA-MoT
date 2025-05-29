# Define the input grid
input_grid = [
    [7, 6, 6, 6, 6],
    [6, 7, 1, 6, 6],
    [6, 1, 9, 9, 6],
    [6, 1, 6, 8, 8]
]

# Duplicate the grid by appending it below itself
output_grid = input_grid + input_grid

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))