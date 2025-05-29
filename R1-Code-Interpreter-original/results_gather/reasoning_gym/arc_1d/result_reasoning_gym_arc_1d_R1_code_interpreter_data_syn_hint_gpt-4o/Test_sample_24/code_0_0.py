def transform_grid(input_grid):
    # Count leading zeros
    leading_zeros = 0
    for num in input_grid:
        if num == 0:
            leading_zeros += 1
        else:
            break
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    non_zero_start = leading_zeros
    for i in range(non_zero_start, len(input_grid)):
        if input_grid[i] != 0:
            output_grid[i] = input_grid[i]
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 0, 0, 0]
# Transform the test input
test_output = transform_grid(test_input)
print(test_output)