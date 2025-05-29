def transform_grid(input_grid):
    # Find the first non-zero element
    start = 0
    while start < len(input_grid) and input_grid[start] == 0:
        start += 1
    
    # Find the last non-zero element
    end = len(input_grid) - 1
    while end >= 0 and input_grid[end] == 0:
        end -= 1
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[start:end+1]
    
    # Calculate the number of leading and trailing zeros
    leading_zeros = start
    trailing_zeros = len(input_grid) - end - 1
    
    # Create the output grid
    output_grid = [0] * leading_zeros + non_zero_sequence + [0] * trailing_zeros
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
output = transform_grid(test_input)
print(output)