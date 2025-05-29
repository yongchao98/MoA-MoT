def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If there is no non-zero number, return the input grid as is
    if first_non_zero_index is None:
        return input_grid
    
    # Extract the sequence of non-zero numbers
    non_zero_sequence = input_grid[first_non_zero_index:]
    
    # Calculate the number of leading zeros in the output
    leading_zeros = first_non_zero_index - 1
    
    # Construct the output grid
    output_grid = [0] * leading_zeros + non_zero_sequence + [0] * (len(input_grid) - len(non_zero_sequence) - leading_zeros)
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)