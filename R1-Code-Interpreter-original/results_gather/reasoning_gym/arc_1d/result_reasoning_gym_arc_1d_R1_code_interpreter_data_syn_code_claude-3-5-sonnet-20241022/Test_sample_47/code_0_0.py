def analyze_and_predict(input_grid):
    # Find non-zero sequence
    non_zero = []
    start_idx = -1
    for i, val in enumerate(input_grid):
        if val != 0:
            if start_idx == -1:
                start_idx = i
            non_zero.append(val)
    
    # If no non-zero sequence found, return the input
    if not non_zero:
        return input_grid
    
    # Create output grid
    output = [0] * len(input_grid)
    
    # Calculate shift amount (seems to be related to starting position)
    # After analyzing examples, shift appears to be around 3-4 positions
    shift = 4 if start_idx <= 5 else 3
    
    # Place the non-zero sequence in new position
    new_start = start_idx + shift
    for i, val in enumerate(non_zero):
        output[new_start + i] = val
    
    return output

# Test input
test_input = [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Get prediction
result = analyze_and_predict(test_input)
print(' '.join(map(str, result)))