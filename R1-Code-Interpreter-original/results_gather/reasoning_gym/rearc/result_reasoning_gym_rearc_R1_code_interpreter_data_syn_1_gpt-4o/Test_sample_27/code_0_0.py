def transform_grid(input_grid):
    # Find unique non-zero numbers and their positions
    positions = {}
    for r, row in enumerate(input_grid):
        for c, value in enumerate(row):
            if value != 0:
                if value not in positions:
                    positions[value] = {'min_row': r, 'max_row': r, 'min_col': c, 'max_col': c}
                else:
                    positions[value]['min_row'] = min(positions[value]['min_row'], r)
                    positions[value]['max_row'] = max(positions[value]['max_row'], r)
                    positions[value]['min_col'] = min(positions[value]['min_col'], c)
                    positions[value]['max_col'] = max(positions[value]['max_col'], c)
    
    # Create the output grid
    output_grid = [[0 for _ in range(len(input_grid[0]))] for _ in range(len(input_grid))]
    
    # Fill the output grid based on the positions
    for value, pos in positions.items():
        for r in range(pos['min_row'], pos['max_row'] + 1):
            for c in range(pos['min_col'], pos['max_col'] + 1):
                output_grid[r][c] = value
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 3, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))