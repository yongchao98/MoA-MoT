def apply_rule(input_grid):
    # Determine the size of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0]) if rows > 0 else 0
    
    # Create the output grid with an additional border
    output_grid = []
    
    # Add the top border
    output_grid.append([0] * (cols + 2))
    
    # Add the left and right borders to each row of the input grid
    for row in input_grid:
        output_grid.append([0] + row + [0])
    
    # Add the bottom border
    output_grid.append([0] * (cols + 2))
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0],
    [0, 0]
]

# Apply the rule
output_grid = apply_rule(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))