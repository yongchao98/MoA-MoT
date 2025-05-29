def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is None:
        # If there are no non-zero elements, return the input as is
        return input_grid
    
    # Shift the sequence of non-zero elements to the left by one position
    output_grid = input_grid[first_non_zero_index:first_non_zero_index + input_grid.count(input_grid[first_non_zero_index])]
    output_grid = [0] * first_non_zero_index + output_grid[1:] + [0]
    
    return output_grid

# Test input
test_input = [8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)