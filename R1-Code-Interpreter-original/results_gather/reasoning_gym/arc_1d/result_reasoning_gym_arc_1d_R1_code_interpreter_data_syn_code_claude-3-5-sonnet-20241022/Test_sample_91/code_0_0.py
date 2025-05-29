def shift_pattern(input_grid):
    # Convert input string to list of numbers
    grid = [int(x) for x in input_grid.split()]
    
    # Find non-zero sequence
    non_zeros = []
    for num in grid:
        if num != 0:
            non_zeros.append(num)
    
    # Create output grid (all zeros)
    output = [0] * 18
    
    # Place non-zero sequence shifted 3 positions right
    start_pos = 0
    for i in range(len(grid)):
        if grid[i] != 0:
            start_pos = i
            break
    
    # Place numbers shifted 3 positions to the right
    new_start = start_pos + 3
    for i in range(len(non_zeros)):
        if new_start + i < len(output):
            output[new_start + i] = non_zeros[i]
    
    # Convert to string
    return ' '.join(map(str, output))

# Test input
test_input = "7 2 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(test_input)
print(shift_pattern(test_input))