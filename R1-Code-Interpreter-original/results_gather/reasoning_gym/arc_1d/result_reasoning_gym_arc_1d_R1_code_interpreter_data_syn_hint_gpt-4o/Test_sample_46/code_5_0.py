def transform_grid(input_grid):
    # Find the last non-zero element
    last_non_zero = input_grid[-1]
    
    # Shift elements left by one position
    transformed_grid = input_grid[1:] + [last_non_zero]
    
    return transformed_grid

# Test input
test_input = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# Apply the transformation
output_grid = transform_grid(test_input)
print(output_grid)