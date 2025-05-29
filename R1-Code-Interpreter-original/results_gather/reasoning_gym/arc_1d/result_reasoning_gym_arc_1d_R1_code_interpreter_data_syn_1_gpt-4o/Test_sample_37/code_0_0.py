def transform_grid(input_grid):
    # Find the first non-zero index
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    # Find the last non-zero index
    last_non_zero = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero:last_non_zero + 1]
    
    # Create the output grid with zeros
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence one position to the right
    output_grid[first_non_zero + 1:first_non_zero + 1 + len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
test_input = [0, 0, 4, 3, 5, 4, 9, 0, 0, 0, 0, 0, 0]
# Transform the test input
test_output = transform_grid(test_input)
print(test_output)