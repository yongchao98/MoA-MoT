def transform_grid(input_grid):
    # Find the first non-zero sequence
    non_zero_sequence = []
    for num in input_grid:
        if num != 0:
            non_zero_sequence.append(num)
        elif non_zero_sequence:
            break
    
    # Remove the non-zero sequence from the input grid
    remaining_grid = [num for num in input_grid if num not in non_zero_sequence]
    
    # Construct the output grid
    output_grid = remaining_grid + non_zero_sequence
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 2, 7, 0]
output = transform_grid(test_input)
print(output)