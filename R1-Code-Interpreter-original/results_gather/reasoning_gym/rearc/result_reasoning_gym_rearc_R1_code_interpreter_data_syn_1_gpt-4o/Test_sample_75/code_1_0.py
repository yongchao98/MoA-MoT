def transform_grid(input_grid):
    n = len(input_grid)
    m = len(input_grid[0])
    
    # Determine the border and inner values
    border_value = input_grid[0][0]
    inner_value = None
    
    # Find the first different value from the border
    for i in range(1, n-1):
        for j in range(1, m-1):
            if input_grid[i][j] != border_value:
                inner_value = input_grid[i][j]
                break
        if inner_value is not None:
            break
    
    # If the grid is uniform, return it as is
    if inner_value is None:
        return input_grid
    
    # Create the output grid
    output_grid = [[inner_value] * m for _ in range(n)]
    
    # Set the border to the original inner value
    for i in range(n):
        output_grid[i][0] = border_value
        output_grid[i][m-1] = border_value
    for j in range(m):
        output_grid[0][j] = border_value
        output_grid[n-1][j] = border_value
    
    # Preserve the center area
    for i in range(1, n-1):
        for j in range(1, m-1):
            if input_grid[i][j] != inner_value:
                output_grid[i][j] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 1, 1, 1, 1, 6, 6],
    [6, 6, 1, 2, 2, 1, 6, 6],
    [6, 6, 1, 2, 2, 1, 6, 6],
    [6, 6, 1, 1, 1, 1, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))