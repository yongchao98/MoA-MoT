def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is None:
        # If there are no non-zero elements, return the input as is
        return input_grid
    
    # Find the last non-zero element
    last_non_zero_index = len(input_grid) - 1 - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)
    
    # Extract the sequence of non-zero elements
    non_zero_sequence = input_grid[first_non_zero_index:last_non_zero_index + 1]
    
    # Shift the sequence of non-zero elements to the left by one position
    shifted_sequence = non_zero_sequence[1:] + [0]
    
    # Construct the output grid
    output_grid = input_grid[:first_non_zero_index] + shifted_sequence + input_grid[last_non_zero_index + 1:]
    
    return output_grid

# Test input
test_input = [8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)