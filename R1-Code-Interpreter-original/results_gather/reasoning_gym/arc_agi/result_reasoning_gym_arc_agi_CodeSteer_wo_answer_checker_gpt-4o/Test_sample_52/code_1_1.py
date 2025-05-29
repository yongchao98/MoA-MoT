def transform_grid(input_grid):
    # Get the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    
    # Iterate over each cell in the grid
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 2:
                # Directly adjacent cells
                for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        output_grid[nr][nc] = 4
                
                # Extended influence (diagonal and two steps away)
                for dr, dc in [(-1, -1), (-1, 1), (1, -1), (1, 1), (-2, 0), (2, 0), (0, -2), (0, 2)]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        output_grid[nr][nc] = 4
    
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

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))