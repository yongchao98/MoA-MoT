def map_input_to_output(input_grid):
    # Define the positions of the top-left corners of the 3x3 subgrids in the input grid
    subgrid_positions = [(0, 0), (0, 4), (0, 8),
                         (4, 0), (4, 4), (4, 8),
                         (8, 0), (8, 4), (8, 8)]
    
    # Initialize the output grid
    output_grid = []
    
    # Iterate over each subgrid position
    for pos in subgrid_positions:
        x, y = pos
        # Check if there is a 5 in the 3x3 subgrid
        contains_five = any(input_grid[i][j] == 5 for i in range(x, x+3) for j in range(y, y+3))
        # Append 6 if there is a 5, otherwise append 3
        output_grid.append(6 if contains_five else 3)
    
    # Reshape the output grid into a 3x3 matrix
    output_grid = [output_grid[i:i+3] for i in range(0, 9, 3)]
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 3, 2, 3, 5, 3, 2, 3, 3, 3],
    [5, 3, 3, 2, 3, 3, 3, 2, 3, 5, 3],
    [3, 3, 5, 2, 3, 3, 3, 2, 3, 3, 3],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [3, 3, 3, 2, 5, 3, 3, 2, 3, 3, 5],
    [3, 3, 5, 2, 3, 5, 3, 2, 3, 3, 3],
    [3, 3, 3, 2, 3, 3, 3, 2, 3, 5, 3],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3],
    [3, 5, 3, 2, 3, 3, 5, 2, 3, 3, 3],
    [3, 3, 3, 2, 3, 3, 3, 2, 5, 3, 5]
]

# Get the output grid
output_grid = map_input_to_output(input_grid)

# Print the output grid
print("<<<", output_grid, ">>>")