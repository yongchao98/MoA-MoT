def transform_grid(input_grid):
    # Determine the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Initialize the output grid
    output_grid = [[0] * num_cols for _ in range(num_rows)]
    
    # Transform each column based on unique values
    for col in range(num_cols):
        # Get unique values in the column
        unique_values = sorted(set(row[col] for row in input_grid))
        # Create a mapping from unique values to a sequence starting from 1
        value_to_sequence = {value: i + 1 for i, value in enumerate(unique_values)}
        # Apply the transformation to the column
        for row in range(num_rows):
            output_grid[row][col] = value_to_sequence[input_grid[row][col]]
    
    return output_grid

# Test input grid
input_grid = [
    [5, 8, 5, 5, 5],
    [5, 8, 5, 5, 5],
    [5, 8, 5, 5, 9],
    [5, 8, 7, 5, 9],
    [1, 8, 7, 5, 9]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))