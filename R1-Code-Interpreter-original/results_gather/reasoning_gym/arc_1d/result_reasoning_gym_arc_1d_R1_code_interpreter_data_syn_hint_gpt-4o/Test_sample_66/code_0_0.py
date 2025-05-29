def transform_grid(input_grid):
    # Find the first non-zero element from the right
    for i in range(len(input_grid) - 1, -1, -1):
        if input_grid[i] != 0:
            break
    
    # If the sequence is at the end, shift it left
    if i == len(input_grid) - 1:
        return input_grid
    
    # Otherwise, shift the sequence left and replace trailing with 0s
    non_zero_sequence = input_grid[i:]
    transformed_grid = input_grid[:i] + non_zero_sequence + [0] * (len(input_grid) - len(non_zero_sequence) - i)
    return transformed_grid

input_grid = [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)