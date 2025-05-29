def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
    # Find the last non-zero element
    last_non_zero = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), len(input_grid))
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    for i in range(first_non_zero, last_non_zero):
        output_grid[i] = input_grid[i]
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 0, 0, 0, 0, 0]
# Transform the test input
test_output = transform_grid(test_input)
print(test_output)