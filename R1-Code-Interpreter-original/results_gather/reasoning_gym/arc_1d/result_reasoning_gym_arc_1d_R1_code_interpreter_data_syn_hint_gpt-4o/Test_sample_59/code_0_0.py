def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If no non-zero number is found, return the input as is
    if first_non_zero_index is None:
        return input_grid
    
    # Extract the sequence of non-zero numbers
    non_zero_sequence = [x for x in input_grid if x != 0]
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence starting from the first non-zero index
    output_grid[first_non_zero_index:first_non_zero_index + len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# Get the output
test_output = transform_grid(test_input)
print(test_output)