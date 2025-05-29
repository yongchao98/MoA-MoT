def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If no non-zero element is found, return the input as is
    if first_non_zero_index is None:
        return input_grid
    
    # Determine the length of the non-zero sequence
    non_zero_sequence = [x for x in input_grid if x != 0]
    sequence_length = len(non_zero_sequence)
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence shifted to the left
    output_grid[first_non_zero_index - sequence_length:first_non_zero_index] = non_zero_sequence
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0]
output = transform_grid(test_input)
print(output)