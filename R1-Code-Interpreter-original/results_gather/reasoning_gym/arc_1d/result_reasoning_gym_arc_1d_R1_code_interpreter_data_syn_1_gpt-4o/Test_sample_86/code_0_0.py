def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is not None:
        # Move the first non-zero element to the start
        output_grid = [0] * first_non_zero_index + input_grid[first_non_zero_index:]
    else:
        # If all are zeros, return the same grid
        output_grid = input_grid
    
    return output_grid

# Test input
test_input = [8, 1, 2, 0, 0, 0, 0, 0, 0, 0, 8, 7, 2, 8, 4, 5, 1, 8, 2]
output = transform_grid(test_input)
print(output)