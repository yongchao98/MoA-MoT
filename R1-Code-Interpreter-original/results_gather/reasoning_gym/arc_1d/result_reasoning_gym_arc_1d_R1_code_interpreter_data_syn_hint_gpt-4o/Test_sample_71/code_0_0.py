def transform_grid(input_grid):
    # Find the position of '2'
    pos_2 = input_grid.index(2)
    
    # Extract the sequence of non-zero numbers other than '2'
    sequence = [x for x in input_grid if x != 0 and x != 2]
    
    # Create the output grid initialized with zeros
    output_grid = [0] * len(input_grid)
    
    # Place '2' in its original position
    output_grid[pos_2] = 2
    
    # Place the sequence in the output grid, shifted to the left
    for i, num in enumerate(sequence):
        output_grid[i] = num
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 2, 0, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0]
# Get the output
output = transform_grid(test_input)
print(output)