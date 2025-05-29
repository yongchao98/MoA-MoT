def extract_subgrid(input_grid):
    n = len(input_grid)
    m = len(input_grid[0])
    
    # Determine the bottom boundary
    bottom = n
    for i in range(n - 1, -1, -1):
        if 9 in input_grid[i]:
            bottom = i + 1
            break
    
    # Determine the right boundary
    right = m
    for j in range(m - 1, -1, -1):
        if any(input_grid[i][j] == 9 for i in range(n)):
            right = j + 1
            break
    
    # Extract the subgrid
    subgrid = [row[:right] for row in input_grid[:bottom]]
    
    return subgrid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 3, 0, 1, 3, 9],
    [9, 9, 9, 9, 9, 0, 3, 0, 0, 9],
    [9, 0, 0, 3, 9, 1, 0, 1, 3, 9],
    [9, 3, 0, 1, 9, 9, 9, 9, 9, 9],
    [9, 1, 0, 0, 9, 9, 9, 9, 9, 9],
    [9, 0, 3, 1, 9, 9, 0, 0, 0, 0],
    [9, 0, 0, 0, 9, 9, 3, 0, 1, 3],
    [9, 3, 0, 1, 9, 9, 0, 3, 0, 0],
    [9, 9, 9, 9, 9, 9, 0, 0, 0, 1]
]

output_grid = extract_subgrid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))