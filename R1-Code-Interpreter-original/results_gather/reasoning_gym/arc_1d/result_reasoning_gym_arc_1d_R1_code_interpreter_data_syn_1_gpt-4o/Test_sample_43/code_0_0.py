def transform_grid(input_grid):
    # Find the first and last non-zero indices
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    last_non_zero = next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)
    if last_non_zero is not None:
        last_non_zero = len(input_grid) - 1 - last_non_zero
    
    # Reverse the segment between first and last non-zero indices
    if first_non_zero is not None and last_non_zero is not None:
        input_grid[first_non_zero:last_non_zero+1] = input_grid[first_non_zero:last_non_zero+1][::-1]
    
    return input_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)