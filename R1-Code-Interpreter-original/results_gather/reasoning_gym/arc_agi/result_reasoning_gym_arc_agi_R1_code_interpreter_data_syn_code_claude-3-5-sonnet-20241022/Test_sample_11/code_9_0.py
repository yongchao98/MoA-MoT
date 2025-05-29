def process_grid(grid):
    # Convert input grid to integers if they're strings
    grid = [[int(x) for x in row] for row in grid]
    
    # Create new grid
    new_grid = []
    for row in grid:
        new_grid.append(row.copy())
    
    # Process zeros
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            # Find zeros and replace them
            if grid[i][j] == 0:
                # Upper half (including middle) uses first row value
                if i <= 4:
                    new_grid[i][j] = grid[0][j]
                # Lower half uses last row value
                else:
                    new_grid[i][j] = grid[-1][j]
    
    # Print for debugging
    for row in new_grid:
        print(' '.join(str(x) for x in row))
    
    return new_grid

# Input grid
grid = [
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,0,8,8,8,0,8,8,6],
    [4,8,8,8,0,8,8,8,8,6],
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,8,0,8,8,8,8,8,6],
    [4,8,8,8,8,8,0,8,8,6],
    [4,8,8,0,8,8,8,8,8,6],
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,8,8,8,0,8,8,0,6],
    [4,0,8,8,8,8,8,8,8,6]
]

result = process_grid(grid)