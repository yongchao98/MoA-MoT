def transform_grid(input_grid):
    output = input_grid.copy()
    current_num = None
    start_idx = None
    
    # Process the grid
    for i in range(len(input_grid)):
        if input_grid[i] != 0:
            if current_num is None:  # Start of new sequence
                current_num = input_grid[i]
                start_idx = i
            elif input_grid[i] != current_num:  # End of sequence, start of new one
                if i - start_idx > 1:  # If sequence was longer than 1
                    output[start_idx+1:i] = [0] * (i - start_idx - 1)
                current_num = input_grid[i]
                start_idx = i
        else:  # Hit a zero
            if current_num is not None:  # End of sequence
                if i - start_idx > 1:  # If sequence was longer than 1
                    output[start_idx+1:i-1] = [0] * (i - start_idx - 2)
                current_num = None
    
    # Handle case where sequence ends at the end of the grid
    if current_num is not None and len(input_grid) - start_idx > 1:
        output[start_idx+1:len(input_grid)-1] = [0] * (len(input_grid) - start_idx - 2)
    
    # Print the result
    print(' '.join(map(str, output)))

# Test input
test_input = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0]
transform_grid(test_input)