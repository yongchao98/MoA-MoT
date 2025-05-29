# Define the input grid
input_grid = [
    [1, 1, 8, 8, 8, 1],
    [1, 1, 8, 7, 8, 1],
    [1, 1, 8, 8, 8, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1]
]

# Define the output grid based on the corrected rule
output_grid = [
    [8, 8, 8],
    [8, 7, 8],
    [8, 8, 8]
]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))