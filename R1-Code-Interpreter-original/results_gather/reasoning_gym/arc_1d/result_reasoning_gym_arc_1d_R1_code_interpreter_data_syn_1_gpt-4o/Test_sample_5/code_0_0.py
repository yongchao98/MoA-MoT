def transform_grid(input_grid):
    # Find the first non-zero index
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero_index:]
    
    # Count the number of zeros at the end of the non-zero sequence
    trailing_zeros_count = len(input_grid) - len(non_zero_sequence)
    
    # Create the output grid
    output_grid = non_zero_sequence + [0] * trailing_zeros_count
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Transform the test input
output = transform_grid(test_input)

# Print the output
print(output)