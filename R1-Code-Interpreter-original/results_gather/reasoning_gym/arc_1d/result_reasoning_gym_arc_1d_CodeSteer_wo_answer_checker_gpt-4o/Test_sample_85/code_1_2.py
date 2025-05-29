def transform_grid(input_grid):
    # Find the first non-zero sequence
    non_zero_sequence = []
    for num in input_grid:
        if num != 0:
            non_zero_sequence.append(num)
        elif non_zero_sequence:
            break
    
    # Create the output grid
    output_grid = non_zero_sequence + [0] * (len(input_grid) - len(non_zero_sequence))
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output
print("<<<" + " ".join(map(str, output_grid)) + ">>>")