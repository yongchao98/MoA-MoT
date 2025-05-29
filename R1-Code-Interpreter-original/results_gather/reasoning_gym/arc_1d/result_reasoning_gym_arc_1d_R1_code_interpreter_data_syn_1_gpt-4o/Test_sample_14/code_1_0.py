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
    
    # Calculate the number of zeros to maintain
    total_zeros = len(input_grid) - len(non_zero_sequence)
    
    # Create the output grid with the non-zero sequence at the leftmost position
    output_grid = non_zero_sequence + [0] * total_zeros
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
output = transform_grid(test_input)
print(output)