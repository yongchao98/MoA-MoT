def transform_grid(input_grid):
    # Define the transformation rules
    transformation_rules = {8: 9, 9: 8, 3: 4, 4: 3, 0: 8, 8: 0}
    
    # Apply the transformation rules to the input grid
    output_grid = []
    for row in input_grid:
        transformed_row = [transformation_rules.get(value, value) for value in row]
        output_grid.append(transformed_row)
    
    return output_grid

# Test input grid
input_grid = [
    [0, 1, 1, 2, 2, 2, 2, 2, 2, 2],
    [1, 0, 1, 2, 1, 1, 1, 2, 2, 2],
    [1, 1, 1, 2, 2, 2, 2, 2, 2, 2],
    [1, 5, 1, 1, 1, 1, 1, 1, 1, 1],
    [5, 5, 5, 1, 1, 0, 1, 0, 1, 1],
    [1, 1, 2, 2, 2, 2, 2, 1, 1, 0],
    [1, 1, 2, 2, 2, 2, 2, 1, 1, 1],
    [1, 1, 2, 2, 1, 2, 2, 6, 6, 6],
    [0, 1, 2, 1, 1, 1, 2, 1, 1, 1],
    [1, 1, 2, 2, 2, 2, 2, 1, 1, 0]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))