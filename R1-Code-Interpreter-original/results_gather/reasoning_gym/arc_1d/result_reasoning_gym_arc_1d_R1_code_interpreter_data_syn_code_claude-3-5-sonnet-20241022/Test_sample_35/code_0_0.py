def transform_grid(input_grid):
    # Find non-zero sequence
    non_zeros = [x for x in input_grid if x != 0]
    
    # Create output grid
    output = []
    
    # Add leading zero
    output.append(0)
    
    # Add non-zero sequence
    output.extend(non_zeros)
    
    # Fill remaining positions with zeros
    while len(output) < len(input_grid):
        output.append(0)
    
    return output

# Test input
test_input = [0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
result = transform_grid(test_input)
print(' '.join(map(str, result)))