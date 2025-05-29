# Define the input grid
input_grid = [
    [6, 3, 0, 6, 3, 9, 9, 9, 9, 9, 9, 0],
    [3, 2, 0, 3, 2, 9, 9, 9, 9, 9, 9, 0],
    [0, 0, 0, 0, 0, 9, 9, 9, 9, 9, 9, 0],
    [9, 9, 9, 6, 3, 9, 9, 9, 9, 9, 9, 0],
    [9, 9, 9, 3, 2, 9, 9, 9, 9, 9, 9, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 3, 0, 6, 3, 0, 6, 3, 0, 6, 3, 0],
    [3, 2, 0, 3, 2, 0, 3, 2, 0, 3, 2, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 3, 0, 6, 3, 0, 6, 3, 0, 6, 3, 0],
    [3, 2, 0, 3, 2, 0, 3, 2, 0, 3, 2, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Initialize the output grid
output_grid = []

# First three rows remain unchanged
output_grid.extend(input_grid[:3])

# Middle rows (4th to 9th) are transformed
pattern = input_grid[:3]
for i in range(3, 9):
    new_row = pattern[i % 3]
    output_grid.append(new_row)

# Last three rows remain unchanged
output_grid.extend(input_grid[9:])

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))