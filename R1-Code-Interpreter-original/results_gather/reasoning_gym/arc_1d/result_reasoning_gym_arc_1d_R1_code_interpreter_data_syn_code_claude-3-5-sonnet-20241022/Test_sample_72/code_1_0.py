def transform_grid(input_grid):
    # Initialize output with zeros
    output = [0] * 16
    
    # Find the non-zero sequence
    non_zero = [(i, val) for i, val in enumerate(input_grid) if val != 0]
    
    # Apply the transformation rule:
    # Shift each non-zero value 3 positions to the right
    for i, (pos, val) in enumerate(non_zero):
        new_pos = pos + 3
        if new_pos < len(output):
            output[new_pos] = val
    
    return output

# Test input
test_input = [0,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0]
result = transform_grid(test_input)
print(' '.join(map(str, result)))