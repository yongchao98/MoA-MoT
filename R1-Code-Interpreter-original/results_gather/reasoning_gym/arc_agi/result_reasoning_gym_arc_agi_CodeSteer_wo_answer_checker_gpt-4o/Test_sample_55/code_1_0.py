def transform_grid(input_grid):
    """
    Transforms the input grid by reversing the order of the rows.
    This function can be extended to apply more complex transformations if needed.
    """
    # Reverse the order of the rows
    return input_grid[::-1]

# Define the input grid
input_grid = [
    [4, 6, 6],
    [6, 4, 3],
    [0, 6, 3]
]

# Apply the transformation
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))