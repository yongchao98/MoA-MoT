def transform_grid(input_grid):
    # Find the first and last non-zero elements
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    last_non_zero_index = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1
    
    # If no non-zero element is found, return the input as is
    if first_non_zero_index is None or last_non_zero_index is None:
        return input_grid
    
    # Determine the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero_index:last_non_zero_index + 1]
    
    # Calculate the starting position to center the non-zero sequence
    total_length = len(input_grid)
    non_zero_length = len(non_zero_sequence)
    start_position = (total_length - non_zero_length) // 2
    
    # Create the output grid
    output_grid = [0] * total_length
    
    # Place the non-zero sequence in the calculated position
    for i, value in enumerate(non_zero_sequence):
        output_grid[start_position + i] = value
    
    return output_grid

# Test input
test_input = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# Transform the test input
test_output = transform_grid(test_input)
print(test_output)