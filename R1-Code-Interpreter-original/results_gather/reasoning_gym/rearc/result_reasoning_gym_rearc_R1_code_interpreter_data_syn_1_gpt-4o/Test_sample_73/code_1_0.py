def transform_grid(input_grid):
    # Determine the columns to retain based on the presence of specific values
    retained_columns = []
    for col in range(len(input_grid[0])):
        if any(input_grid[row][col] in [1, 4, 6, 8] for row in range(len(input_grid))):
            retained_columns.append(col)
    
    # Create the output grid based on the retained columns and transformation rule
    output_grid = []
    for row in input_grid:
        new_row = []
        for col in retained_columns:
            if row[col] in [1, 4, 6, 8]:
                new_row.append(3)
            else:
                new_row.append(row[col])
        output_grid.append(new_row)
    
    return output_grid

# Test input grid
test_input_grid = [
    [9, 7, 8, 7, 7],
    [7, 7, 8, 7, 1]
]

# Transform the test input grid
predicted_output_grid = transform_grid(test_input_grid)

# Print the predicted output grid
for row in predicted_output_grid:
    print(" ".join(map(str, row)))