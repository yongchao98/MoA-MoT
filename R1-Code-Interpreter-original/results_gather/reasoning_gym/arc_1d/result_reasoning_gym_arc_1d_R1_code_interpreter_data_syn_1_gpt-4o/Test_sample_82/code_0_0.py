def transform_grid(input_grid):
    # Find the first non-zero sequence
    non_zero_sequence = []
    for num in input_grid:
        if num != 0:
            non_zero_sequence.append(num)
        elif non_zero_sequence:
            break
    
    # Calculate the number of leading zeros
    leading_zeros = input_grid.index(non_zero_sequence[0])
    
    # Move the sequence to the left by the number of leading zeros
    output_grid = [0] * len(input_grid)
    output_grid[leading_zeros:leading_zeros+len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
test_input = [6, 8, 4, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 2]
output = transform_grid(test_input)
print(output)