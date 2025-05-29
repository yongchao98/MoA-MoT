# Define the input grid
input_grid = [
    [7, 0, 0],
    [7, 7, 9],
    [7, 7, 7],
    [7, 5, 9],
    [7, 9, 7],
    [7, 7, 7],
    [9, 6, 6],
    [7, 7, 9],
    [7, 7, 7],
    [7, 4, 9],
    [7, 4, 7]
]

# Apply the transformation rule
output_grid = [[row[1], row[0], row[2]] for row in input_grid]

# Print the output grid
print(output_grid)