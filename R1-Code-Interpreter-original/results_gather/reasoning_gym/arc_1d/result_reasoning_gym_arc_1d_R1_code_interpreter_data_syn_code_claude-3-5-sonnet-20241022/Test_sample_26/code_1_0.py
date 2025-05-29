def process_grid(grid):
    # Find the first sequence of zeros
    start_zero = -1
    end_zero = -1
    in_zero_sequence = False
    
    for i in range(len(grid)):
        if grid[i] == 0 and not in_zero_sequence:
            start_zero = i
            in_zero_sequence = True
        elif grid[i] != 0 and in_zero_sequence:
            end_zero = i - 1
            break
        if in_zero_sequence and i == len(grid) - 1:
            end_zero = i
    
    if start_zero == -1:
        return grid
    
    # Create output grid
    result = grid.copy()
    zero_length = end_zero - start_zero + 1
    new_start = max(0, start_zero - 2)  # Shift left by 2
    
    # Reset the original zero positions
    for i in range(start_zero, end_zero + 1):
        if i + zero_length < len(grid):
            result[i] = grid[i + zero_length]
        else:
            result[i] = grid[i]
    
    # Place zeros in new position
    for i in range(new_start, new_start + zero_length):
        result[i] = 0
    
    print(' '.join(map(str, result)))

# Test input
test_input = [2, 2, 2, 2, 2, 2, 2, 0, 0, 2, 2, 2, 2, 2, 2, 2]
process_grid(test_input)