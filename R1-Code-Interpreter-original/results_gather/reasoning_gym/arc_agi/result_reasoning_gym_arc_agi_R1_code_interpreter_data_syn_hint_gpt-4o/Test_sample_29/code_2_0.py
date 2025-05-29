def flood_fill(grid, x, y, target, replacement):
    if x < 0 or x >= len(grid) or y < 0 or y >= len(grid[0]):
        return
    if grid[x][y] != target:
        return
    grid[x][y] = replacement
    flood_fill(grid, x + 1, y, target, replacement)
    flood_fill(grid, x - 1, y, target, replacement)
    flood_fill(grid, x, y + 1, target, replacement)
    flood_fill(grid, x, y - 1, target, replacement)

def transform_grid(input_grid):
    # Initialize the output grid with the same dimensions as the input grid
    output_grid = [[3] * len(row) for row in input_grid]
    
    # Iterate over each cell in the grid
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            value = input_grid[i][j]
            if value in [9, 4, 6, 5, 0, 3, 7]:
                output_grid[i][j] = value
            elif value == 2:
                # Use flood fill to transform contiguous blocks of 2 into 1
                flood_fill(output_grid, i, j, 3, 1)
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 8, 8],
    [8, 8, 8, 0, 8, 8, 8, 2, 2, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 8],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8],
    [8, 8, 8, 8, 8, 2, 8, 8, 8, 8, 8, 2, 8],
    [8, 8, 8, 8, 8, 2, 2, 2, 8, 8, 8, 6, 8],
    [8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 2, 8],
    [8, 8, 8, 8, 8, 8, 8, 2, 8, 4, 8, 2, 8],
    [8, 8, 9, 8, 8, 8, 8, 2, 8, 8, 8, 2, 8],
    [2, 2, 2, 8, 8, 8, 2, 2, 8, 8, 2, 2, 8],
    [8, 8, 2, 2, 8, 8, 2, 8, 8, 8, 2, 8, 8],
    [8, 8, 8, 2, 8, 8, 2, 2, 2, 2, 2, 8, 8],
    [8, 7, 8, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))