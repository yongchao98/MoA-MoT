def transform_grid(input_grid):
    # Get the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    
    # Directions for surrounding cells (8 directions)
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), 
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    # Iterate over each cell in the grid
    for r in range(rows):
        for c in range(cols):
            # If the current cell is zero
            if input_grid[r][c] == 0:
                # Replace surrounding cells with 5
                for dr, dc in directions:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        if input_grid[nr][nc] != 0:
                            output_grid[nr][nc] = 5
    
    return output_grid

# Test input grid
input_grid = [
    [1, 7, 1, 1, 1, 1, 1, 1, 7, 7, 1, 1, 7, 7, 7, 1, 1, 1],
    [1, 1, 1, 1, 1, 7, 1, 1, 7, 1, 1, 1, 7, 0, 0, 7, 0, 7],
    [1, 7, 1, 1, 7, 1, 1, 7, 7, 7, 1, 1, 0, 7, 7, 0, 0, 1],
    [1, 7, 1, 7, 7, 7, 7, 7, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [7, 7, 7, 7, 1, 1, 7, 1, 7, 7, 1, 7, 1, 1, 1, 7, 7, 1],
    [1, 1, 1, 1, 7, 1, 1, 7, 0, 0, 7, 0, 7, 1, 7, 7, 7, 7],
    [1, 1, 1, 1, 7, 7, 1, 7, 0, 0, 0, 0, 7, 1, 1, 1, 1, 1],
    [7, 1, 7, 7, 7, 1, 7, 7, 1, 7, 1, 1, 7, 1, 7, 7, 1, 7],
    [1, 1, 1, 1, 7, 1, 7, 7, 7, 7, 7, 1, 1, 1, 7, 1, 1, 1],
    [7, 1, 1, 1, 1, 7, 7, 1, 1, 7, 7, 1, 7, 1, 7, 1, 7, 7],
    [7, 7, 7, 7, 1, 1, 7, 1, 1, 7, 0, 7, 1, 7, 1, 7, 7, 1],
    [7, 1, 7, 7, 7, 1, 1, 1, 7, 7, 7, 0, 7, 1, 1, 7, 7, 1],
    [1, 7, 7, 7, 1, 1, 1, 7, 1, 1, 7, 7, 7, 1, 7, 7, 1, 7],
    [1, 7, 0, 0, 7, 0, 0, 7, 1, 1, 1, 1, 7, 7, 7, 7, 7, 7],
    [7, 0, 7, 7, 0, 7, 0, 1, 1, 7, 1, 7, 7, 7, 7, 7, 7, 1],
    [1, 7, 0, 0, 7, 0, 0, 1, 7, 1, 7, 7, 1, 1, 7, 1, 7, 1],
    [1, 7, 7, 7, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))