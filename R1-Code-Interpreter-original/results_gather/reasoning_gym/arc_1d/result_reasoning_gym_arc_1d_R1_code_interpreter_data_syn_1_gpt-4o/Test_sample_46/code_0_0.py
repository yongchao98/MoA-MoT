def transform_grid(input_grid):
    # Count trailing zeros
    trailing_zeros = 0
    for num in reversed(input_grid):
        if num == 0:
            trailing_zeros += 1
        else:
            break
    
    # Shift the first `trailing_zeros` elements to the end
    if trailing_zeros > 0:
        return input_grid[trailing_zeros:] + input_grid[:trailing_zeros]
    else:
        return input_grid

# Test input
test_input = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# Transform the test input
test_output = transform_grid(test_input)
print(test_output)