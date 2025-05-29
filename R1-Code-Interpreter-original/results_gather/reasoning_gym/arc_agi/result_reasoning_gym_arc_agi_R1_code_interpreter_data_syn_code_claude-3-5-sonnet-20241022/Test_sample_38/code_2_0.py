def print_grid(grid):
    for row in grid:
        print(' '.join(map(str, row)))

def find_pattern_center(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    for i in range(1, rows-1):
        for j in range(1, cols-1):
            if (grid[i][j] == 6 and
                all(grid[i+di][j+dj] == 9 for di, dj in [(-1,0), (1,0), (0,-1), (0,1)])):
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
        
        # Place the plus (+) pattern in new position
        output_grid[new_i][new_j] = 6  # center
        output_grid[new_i-1][new_j] = 9  # top
        output_grid[new_i+1][new_j] = 9  # bottom
        output_grid[new_i][new_j-1] = 9  # left
        output_grid[new_i][new_j+1] = 9  # right
        
        # Remove the original 6 from where the pattern was
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