def expand_grid(input_grid):
    # Determine the size of the input grid
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    
    # Create the output grid by expanding each element into a 3x3 block
    output_grid = []
    for row in input_grid:
        expanded_rows = [[element] * 3 for element in row]
        for _ in range(3):
            output_grid.append([element for expanded_row in expanded_rows for element in expanded_row])
    
    return output_grid

# Test input grid
input_grid = [
    [2, 2],
    [2, 2]
]

# Generate the output grid
output_grid = expand_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))