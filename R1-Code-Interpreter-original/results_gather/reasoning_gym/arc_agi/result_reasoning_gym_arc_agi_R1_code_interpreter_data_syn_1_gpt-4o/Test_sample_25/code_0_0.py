def transform_grid(input_grid):
    # Define the transformation rules based on observations
    transformation_rules = {
        8: 1,
        3: 1,
        0: 1,
        4: 2,
        6: 6,
        9: 9,
        5: 5
    }
    
    # Initialize the output grid
    output_grid = []
    
    # Apply the transformation rules
    for row in input_grid:
        transformed_row = []
        for value in row:
            if value in transformation_rules:
                transformed_row.append(transformation_rules[value])
            else:
                transformed_row.append(value)
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

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))