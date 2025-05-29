def apply_rule(input_grid):
    output_grid = []
    for row in input_grid:
        # Find the number to propagate
        number_to_propagate = row[-1]
        # Create a new row with the number propagated to the right
        new_row = row[:-1] + [number_to_propagate]
        output_grid.append(new_row)
    return output_grid

# Test input grid
input_grid = [
    [4, 1, 1, 0, 1, 1, 1],
    [1, 4, 1, 0, 1, 1, 1],
    [1, 1, 1, 0, 1, 1, 1]
]

# Apply the rule
output_grid = apply_rule(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))