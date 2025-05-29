def transform_grid(input_grid):
    # Get the non-zero number
    num = next(x for x in input_grid if x != 0)
    
    # Count initial sequence
    initial_count = 0
    for x in input_grid:
        if x == 0:
            break
        initial_count += 1
    
    # Count total zeros
    zero_count = input_grid.count(0)
    
    # Create output grid
    output = [0] * 24  # Initialize with zeros
    
    # Place initial sequence + 2
    for i in range(initial_count + 2):
        output[i] = num
    
    # Keep the last position as it was in input
    output[-1] = input_grid[-1]
    
    # Fill remaining positions after zeros with the number
    non_zero_end = 24 - zero_count - (initial_count + 2)
    if non_zero_end > 0:
        start_idx = initial_count + 2 + zero_count
        for i in range(start_idx, 23):  # Stop before last position
            output[i] = num
    
    print(' '.join(map(str, result)))

# Test input
test = [4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4]
result = transform_grid(test)