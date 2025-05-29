def print_grid(grid):
    for row in grid:
        print(' '.join(map(str, row)))

def find_pattern_center(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    for i in range(1, rows-1):
        for j in range(1, cols-1):
            if (grid[i][j] == 6 and
                grid[i-1][j] == 9 and grid[i+1][j] == 9 and
                grid[i][j-1] == 9 and grid[i][j+1] == 9):
                return (i, j)
    return None

def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[0]*cols for _ in range(rows)]
    
    # Copy all 6's first
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 6:
                output_grid[i][j] = 6
    
    # Find the 3x3 pattern center
    center = find_pattern_center(input_grid)
    if center:
        i, j = center
        # Move pattern up by 2 positions
        new_i = i - 2
        new_j = j
        
        # Remove the original 6 at the new pattern position
        output_grid[new_i][new_j] = 6
        # Place the 3x3 pattern in new position
        output_grid[new_i-1][new_j] = 9
        output_grid[new_i+1][new_j] = 9
        output_grid[new_i][new_j-1] = 9
        output_grid[new_i][new_j+1] = 9
        # Remove the 6 from the original pattern position
        output_grid[i][j] = 0
    
    return output_grid

# Test input
input_grid = [
    [0,0,0,0,6,0,0],
    [0,0,0,0,0,0,0],
    [0,0,0,0,6,0,0],
    [0,0,0,0,0,0,0],
    [0,0,0,0,6,0,0],
    [0,0,0,0,0,0,0],
    [0,0,0,0,6,0,0],
    [0,0,0,9,9,9,0],
    [0,0,0,9,6,9,0],
    [0,0,0,9,9,9,0],
    [0,0,0,0,6,0,0],
    [0,0,0,0,0,0,0],
    [0,0,0,0,6,0,0],
    [0,0,0,0,0,0,0],
    [0,0,0,0,6,0,0],
    [0,0,0,0,0,0,0],
    [0,0,0,0,6,0,0]
]

output_grid = transform_grid(input_grid)
print_grid(output_grid)