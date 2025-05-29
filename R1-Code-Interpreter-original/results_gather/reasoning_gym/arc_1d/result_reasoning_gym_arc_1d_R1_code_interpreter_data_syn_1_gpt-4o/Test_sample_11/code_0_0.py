def transform_grid(input_grid):
    # Find the first non-zero index
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If there is no non-zero element, return the input as is
    if first_non_zero_index is None:
        return input_grid
    
    # Create the output grid
    output_grid = input_grid[:]
    
    # Shift the first non-zero element one position to the right
    if first_non_zero_index < len(input_grid) - 1:
        output_grid[first_non_zero_index] = 0
        output_grid[first_non_zero_index + 1] = input_grid[first_non_zero_index]
    
    return output_grid

# Test input
test_input = [0, 9, 2, 1, 2, 8, 6, 6, 9, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# Transform the test input
test_output = transform_grid(test_input)
print(test_output)