def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[8 for _ in range(cols)] for _ in range(rows)]
    
    # Step 1: Identify positions of '2' and mark them in the output grid
    positions_of_two = []
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 2:
                output_grid[r][c] = 2
                positions_of_two.append((r, c))
    
    # Step 2: Place '4' in a cross pattern around each '2'
    for r, c in positions_of_two:
        if r > 0:
            output_grid[r-1][c] = 4
        if r < rows - 1:
            output_grid[r+1][c] = 4
        if c > 0:
            output_grid[r][c-1] = 4
        if c < cols - 1:
            output_grid[r][c+1] = 4
    
    # Step 3: Fill entire rows and columns with '4' if they contain a '2'
    for r, c in positions_of_two:
        for i in range(cols):
            if output_grid[r][i] != 2:
                output_grid[r][i] = 4
        for i in range(rows):
            if output_grid[i][c] != 2:
                output_grid[i][c] = 4
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 8],
    [8, 2, 8, 8, 8, 8, 8, 8, 8, 8, 2],
    [2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 2, 8, 8, 4, 8, 2, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 2, 8, 8, 8, 8, 2, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 2]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))