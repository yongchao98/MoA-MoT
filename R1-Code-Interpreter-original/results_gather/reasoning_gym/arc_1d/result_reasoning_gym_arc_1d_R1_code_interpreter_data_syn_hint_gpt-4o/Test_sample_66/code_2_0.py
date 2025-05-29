def transform_grid(input_grid):
    # Find the last non-zero element
    last_non_zero_index = -1
    for i in range(len(input_grid) - 1, -1, -1):
        if input_grid[i] != 0:
            last_non_zero_index = i
            break
    
    # If the last non-zero element is not at the end, no transformation is needed
    if last_non_zero_index == len(input_grid) - 1:
        return input_grid
    
    # Find the first non-zero element from the right
    first_non_zero_index = -1
    for i in range(last_non_zero_index, -1, -1):
        if input_grid[i] != 0:
            first_non_zero_index = i
        else:
            break
    
    # Shift the sequence left and replace trailing with 0s
    non_zero_sequence = input_grid[first_non_zero_index:last_non_zero_index + 1]
    transformed_grid = input_grid[:first_non_zero_index] + non_zero_sequence + [0] * (len(input_grid) - len(non_zero_sequence) - first_non_zero_index)
    return transformed_grid

input_grid = [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)