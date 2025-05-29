def transform_grid(input_grid):
    output_grid = input_grid[:]
    n = len(input_grid)
    
    # Find the non-zero number (other than 1) to expand
    for i in range(n):
        if input_grid[i] != 0 and input_grid[i] != 1:
            num_to_expand = input_grid[i]
            break
    
    # Expand the number to fill surrounding zeros
    for i in range(n):
        if input_grid[i] == 0:
            if (i > 0 and input_grid[i-1] == num_to_expand) or (i < n-1 and input_grid[i+1] == num_to_expand):
                output_grid[i] = num_to_expand
    
    return output_grid

# Test input
test_input = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
test_output = transform_grid(test_input)
print(test_output)