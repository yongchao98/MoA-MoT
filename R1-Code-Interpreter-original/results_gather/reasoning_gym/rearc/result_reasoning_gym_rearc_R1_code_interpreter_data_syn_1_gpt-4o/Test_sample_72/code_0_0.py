def apply_rule(input_grid):
    # Determine the repeating pattern
    pattern = input_grid[0]
    
    # Create the output grid by replacing any different blocks with the pattern
    output_grid = []
    for row in input_grid:
        if row == pattern:
            output_grid.append(row)
        else:
            output_grid.append(pattern)
    
    return output_grid

# Test input grid
input_grid = [
    [3, 0, 3, 0, 3, 0, 3, 0],
    [0, 3, 0, 3, 0, 3, 0, 3],
    [3, 0, 3, 0, 3, 0, 3, 0],
    [0, 3, 0, 3, 0, 3, 5, 5],
    [3, 0, 3, 0, 3, 0, 5, 5],
    [0, 3, 0, 3, 0, 3, 0, 3],
    [3, 0, 3, 0, 3, 0, 3, 0],
    [0, 3, 0, 3, 0, 3, 0, 3]
]

# Apply the rule
output_grid = apply_rule(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))