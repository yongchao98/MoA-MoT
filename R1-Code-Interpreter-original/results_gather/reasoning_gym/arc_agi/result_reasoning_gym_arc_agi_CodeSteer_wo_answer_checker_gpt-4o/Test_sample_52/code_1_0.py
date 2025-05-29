def transform_grid(input_grid):
    # Get the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    
    # Function to apply influence of '2' at a given position
    def apply_influence(r, c):
        # Influence pattern: direct neighbors and extended influence
        influence_pattern = [
            (0, 0), (-1, 0), (1, 0), (0, -1), (0, 1),  # Direct neighbors
            (-1, -1), (-1, 1), (1, -1), (1, 1),       # Diagonals
            (-2, 0), (2, 0), (0, -2), (0, 2)          # Two steps away
        ]
        
        for dr, dc in influence_pattern:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                output_grid[nr][nc] = 4
    
    # Iterate over each cell in the grid
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 2:
                apply_influence(r, c)
    
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