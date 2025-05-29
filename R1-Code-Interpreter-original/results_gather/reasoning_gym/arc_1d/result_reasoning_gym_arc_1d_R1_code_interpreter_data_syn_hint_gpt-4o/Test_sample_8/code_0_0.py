def transform_grid(input_grid):
    # Find the sequence of non-zero numbers
    non_zero_sequence = [x for x in input_grid if x != 0 and x != 2]
    
    # Find the position of '2'
    position_of_two = input_grid.index(2)
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence at the beginning
    output_grid[:len(non_zero_sequence)] = non_zero_sequence
    
    # Place '2' at its original position
    output_grid[position_of_two] = 2
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)