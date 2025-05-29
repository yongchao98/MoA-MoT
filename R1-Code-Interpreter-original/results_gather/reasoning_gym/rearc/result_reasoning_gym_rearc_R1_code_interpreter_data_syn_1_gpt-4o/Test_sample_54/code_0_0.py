def extract_subgrid(input_grid, rows, cols):
    return [row[:cols] for row in input_grid[:rows]]

# Test input grid
input_grid = [
    [4, 4, 4, 4, 8, 8, 4, 4, 4],
    [4, 4, 4, 4, 6, 7, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 2, 2, 2, 2, 4, 4],
    [4, 4, 4, 2, 2, 2, 2, 4, 4],
    [4, 4, 4, 2, 2, 2, 2, 4, 4],
    [4, 4, 4, 2, 2, 2, 2, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4]
]

# Extracting a 3x2 subgrid
output_grid = extract_subgrid(input_grid, 3, 2)
for row in output_grid:
    print(' '.join(map(str, row)))