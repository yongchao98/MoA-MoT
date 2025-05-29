def transform_grid(input_grid):
    # Find the longest contiguous sequence of non-zero numbers
    max_start = max_end = start = 0
    max_length = 0
    in_sequence = False
    
    for i, num in enumerate(input_grid):
        if num != 0:
            if not in_sequence:
                start = i
                in_sequence = True
            if i == len(input_grid) - 1 or input_grid[i + 1] == 0:
                end = i
                length = end - start + 1
                if length > max_length:
                    max_length = length
                    max_start = start
                    max_end = end
                in_sequence = False
    
    # Create the output grid
    output_grid = input_grid[:]
    if max_length > 0:
        # Shift the sequence one position to the right
        output_grid[max_start:max_end + 1] = [0] + input_grid[max_start:max_end]
        output_grid[max_start - 1] = input_grid[max_end]
    
    return output_grid

# Test input
test_input = [6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 6]
output = transform_grid(test_input)
print(output)