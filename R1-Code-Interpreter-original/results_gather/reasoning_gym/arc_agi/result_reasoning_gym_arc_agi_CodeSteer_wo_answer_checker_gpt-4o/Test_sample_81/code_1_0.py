def transform_grid(input_grid):
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    output_grid = []
    
    for row in input_grid:
        # Separate 6s and non-6s
        sixes = [x for x in row if x == 6]
        non_sixes = [x for x in row if x != 6]
        
        # Combine them with non-6s at the end
        new_row = sixes + non_sixes
        output_grid.append(new_row)
    
    return output_grid

# Test input grid
input_grid = [
    [6, 2, 6, 2, 2],
    [8, 6, 6, 8, 6],
    [6, 6, 9, 6, 6],
    [1, 6, 6, 1, 6],
    [3, 6, 6, 6, 6]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))