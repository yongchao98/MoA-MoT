def transform_grid(input_grid):
    # Find the number of leading non-zero elements
    non_zero_count = 0
    for num in input_grid:
        if num != 0:
            non_zero_count += 1
        else:
            break
    
    # Create the output grid
    output_grid = [0] * 24
    
    # Shift the non-zero elements to the middle
    output_grid[4:4+non_zero_count] = input_grid[0:non_zero_count]
    
    # Fill the rest of the grid
    output_grid[4+non_zero_count:] = input_grid[non_zero_count:]
    
    return output_grid

# Test input
test_input = [4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4]
output = transform_grid(test_input)
print(output)