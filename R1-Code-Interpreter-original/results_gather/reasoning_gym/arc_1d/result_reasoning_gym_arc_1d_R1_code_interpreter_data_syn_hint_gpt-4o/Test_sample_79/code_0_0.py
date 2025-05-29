def transform_grid(input_grid):
    n = len(input_grid)
    # Find the first and last non-zero elements
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    last_non_zero = next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)
    if last_non_zero is not None:
        last_non_zero = n - last_non_zero - 1
    
    # Create the output grid
    output_grid = [0] * n
    
    # Place the non-zero elements in the new positions
    if first_non_zero is not None and last_non_zero is not None:
        # Calculate the number of zeros to shift
        num_zeros = n - (last_non_zero - first_non_zero + 1)
        left_zeros = num_zeros // 2
        right_zeros = num_zeros - left_zeros
        
        # Fill the output grid
        output_grid[left_zeros:left_zeros + (last_non_zero - first_non_zero + 1)] = input_grid[first_non_zero:last_non_zero + 1]
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)