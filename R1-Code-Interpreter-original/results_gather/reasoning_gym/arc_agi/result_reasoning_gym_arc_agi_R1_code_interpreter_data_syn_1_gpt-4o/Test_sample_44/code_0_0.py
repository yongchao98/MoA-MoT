def transform_grid(input_grid):
    # Extract rows
    row1, row2, row3 = input_grid
    
    # Find the largest number in the grid
    largest = max(max(row1), max(row2), max(row3))
    
    # Apply transformations based on observations
    # First row: Swap first two elements with the last element
    row1_transformed = [row1[2], row1[0], row1[1]]
    
    # Second row: Replace first element with the largest number
    row2_transformed = [largest, row2[1], row2[2]]
    
    # Third row: Remain unchanged
    row3_transformed = row3
    
    return [row1_transformed, row2_transformed, row3_transformed]

# Test input
input_grid = [
    [1, 8, 8],
    [1, 8, 8],
    [8, 1, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
print(output_grid)