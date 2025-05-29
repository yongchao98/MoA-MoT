def process_grid(grid):
    n = len(grid)
    output = [row[:] for row in grid]
    
    # First pass: identify all positions containing 3
    three_positions = []
    for i in range(n):
        for j in range(n):
            if grid[i][j] == 3:
                three_positions.append((i, j))
    
    # For each position with 3, look at surrounding non-3 values
    for i, j in three_positions:
        # Look at values above and below
        above = grid[i-1][j] if i > 0 else None
        below = grid[i+1][j] if i < n-1 else None
        left = grid[i][j-1] if j > 0 else None
        right = grid[i][j+1] if j < n-1 else None
        
        # If any surrounding value is not 3, use that value
        surrounding = [v for v in [above, below, left, right] if v is not None and v != 3]
        if surrounding:
            output[i][j] = surrounding[0]
    
    # Second pass to ensure symmetry
    for i in range(n):
        for j in range(n):
            # Mirror horizontally
            output[i][n-1-j] = output[i][j]
            # Mirror vertically
            output[n-1-i][j] = output[i][j]

    return output

# Parse input grid
input_grid = [
    [4,4,7,5,2,7,9,8,8,9,7,2,5,7,4,4],
    [4,5,5,7,7,7,8,8,8,8,7,7,7,5,5,4],
    [7,5,8,2,9,8,8,3,3,3,8,9,2,8,5,7],
    [5,7,2,5,8,8,5,3,3,3,8,8,5,2,7,5],
    [2,7,3,3,3,3,8,3,3,3,4,7,8,9,7,2],
    [7,7,3,3,3,3,8,0,0,8,4,4,8,8,7,7],
    [9,8,3,3,3,3,0,0,0,0,8,8,5,8,8,9],
    [8,8,3,3,3,3,0,4,4,0,0,8,5,5,8,8],
    [8,8,5,5,8,0,0,4,4,0,0,8,5,5,8,8],
    [9,8,8,5,8,8,0,0,0,0,8,8,5,8,8,9],
    [7,7,8,8,4,4,8,0,0,8,4,4,8,8,7,7],
    [2,7,9,8,7,4,8,8,8,8,4,7,8,9,7,2],
    [5,7,2,5,8,8,5,5,5,5,8,8,5,2,7,5],
    [7,5,8,2,9,8,8,5,5,8,8,9,2,8,5,7],
    [4,5,5,7,7,7,8,8,8,8,7,7,7,5,5,4],
    [4,4,7,5,2,7,9,8,8,9,7,2,5,7,4,4]
]

result = process_grid(input_grid)

# Print the result in the required format
print("<<<")
for row in result:
    print(" ".join(map(str, row)))
print(">>>")