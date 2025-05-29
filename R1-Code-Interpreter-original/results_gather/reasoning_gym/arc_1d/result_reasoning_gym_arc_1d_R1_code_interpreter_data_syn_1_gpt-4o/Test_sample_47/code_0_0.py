def transform_grid(input_grid):
    # Find the first non-zero element
    non_zero_start = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[non_zero_start:]
    
    # Count the number of non-zero elements
    non_zero_count = len(non_zero_sequence)
    
    # Create the output grid
    output_grid = [0] * (len(input_grid) - non_zero_count) + non_zero_sequence
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# Transform the test input
test_output = transform_grid(test_input)
print(test_output)