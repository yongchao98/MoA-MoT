def transform_grid(input_grid):
    output_grid = input_grid[:]
    # Find the first sequence of consecutive '2's
    start_index = None
    for i in range(len(input_grid)):
        if input_grid[i] == 2:
            if start_index is None:
                start_index = i
        else:
            if start_index is not None:
                break
    
    # If a sequence was found, transform it
    if start_index is not None:
        for i in range(start_index, len(input_grid)):
            if input_grid[i] == 2:
                output_grid[i] = 8
            else:
                break
    
    return output_grid

# Test input
test_input = [0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 0]
test_output = transform_grid(test_input)
print(test_output)