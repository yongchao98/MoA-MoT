def transform_grid(input_grid):
    # Get the non-zero number
    num = next(x for x in input_grid if x != 0)
    
    # Count zeros
    zero_count = input_grid.count(0)
    
    # Count initial non-zero sequence
    initial_count = 0
    for x in input_grid:
        if x == 0:
            break
        initial_count += 1
    
    # Create output grid
    output = []
    # Add initial number + 2 more
    output.extend([num] * (initial_count + 2))
    # Add zeros
    output.extend([0] * zero_count)
    # Fill the rest with the number
    output.extend([num] * (24 - len(output)))
    
    return output

# Test input
test = [4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4]
result = transform_grid(test)
print(' '.join(map(str, result)))