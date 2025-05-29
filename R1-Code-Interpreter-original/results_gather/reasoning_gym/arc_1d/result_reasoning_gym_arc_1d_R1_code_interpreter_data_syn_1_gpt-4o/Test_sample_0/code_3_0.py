def transform_grid(input_grid):
    # Find the first and last non-zero indices
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    last_non_zero = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1
    
    # Calculate the new positions for the non-zero sequence
    non_zero_length = last_non_zero - first_non_zero + 1
    new_first = (len(input_grid) - non_zero_length) // 2
    new_last = new_first + non_zero_length - 1
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    output_grid[new_first:new_last+1] = input_grid[first_non_zero:last_non_zero+1]
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 9, 9, 9, 9]
output = transform_grid(test_input)
print(output)