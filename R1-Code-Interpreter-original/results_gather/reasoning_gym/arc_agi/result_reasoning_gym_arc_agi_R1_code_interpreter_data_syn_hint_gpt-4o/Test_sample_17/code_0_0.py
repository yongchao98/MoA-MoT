# Define the input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 8, 8, 4, 4, 4, 8, 8, 4],
    [4, 8, 8, 8, 8, 4, 8, 8, 8, 8],
    [4, 4, 8, 4, 4, 4, 4, 8, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 2, 2, 4],
    [4, 4, 8, 8, 4, 4, 1, 1, 2, 1],
    [4, 8, 8, 8, 8, 4, 4, 1, 4, 4],
    [4, 4, 8, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
]

# Apply the transformation rule
output_grid = []
for row in input_grid:
    new_row = [0 if x == 8 else x for x in row]
    output_grid.append(new_row)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))