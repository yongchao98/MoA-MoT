def find_9_blocks(grid):
    height = len(grid)
    width = len(grid[0])
    blocks = []
    
    # Find all 2x2 blocks of 9's
    for i in range(height-1):
        for j in range(width-1):
            if (grid[i][j] == 9 and grid[i][j+1] == 9 and
                grid[i+1][j] == 9 and grid[i+1][j+1] == 9):
                blocks.append((i, j))
    return blocks

def apply_pattern(grid):
    # Create a copy of the input grid
    output = [row[:] for row in grid]
    
    # Find all 2x2 blocks of 9's
    blocks = find_9_blocks(grid)
    
    # For each block, apply the pattern
    for block_i, block_j in blocks:
        # Check and place 5 (top-left)
        if block_i > 0 and block_j > 0:
            output[block_i-1][block_j-1] = 5
            
        # Check and place 4 (top-right)
        if block_i > 0 and block_j + 2 < len(grid[0]):
            output[block_i-1][block_j+2] = 4
            
        # Check and place 1 (bottom-left)
        if block_i + 2 < len(grid) and block_j > 0:
            output[block_i+2][block_j-1] = 1
            
        # Check and place 8 (bottom-right)
        if block_i + 2 < len(grid) and block_j + 2 < len(grid[0]):
            output[block_i+2][block_j+2] = 8
    
    return output

# Input grid
input_grid = [
    [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3],
    [3,9,9,3,3,3,3,3,3,3,3,3,3,3,3],
    [3,9,9,3,3,3,3,3,3,3,9,9,3,3,3],
    [3,3,3,3,3,3,3,3,3,3,9,9,3,3,3],
    [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3],
    [3,3,3,3,9,9,3,3,3,3,3,3,3,3,3],
    [3,3,3,3,9,9,3,3,3,3,3,3,3,3,3],
    [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3],
    [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3],
    [3,3,3,3,3,3,3,3,3,3,3,9,9,3,3],
    [3,3,3,3,3,3,3,3,3,3,3,9,9,3,3],
    [3,9,9,3,3,3,3,3,3,3,3,3,3,3,3],
    [3,9,9,3,3,3,9,9,3,3,3,3,3,3,3],
    [3,3,3,3,3,3,9,9,3,3,3,3,3,3,3],
    [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]
]

# Get the output grid
output_grid = apply_pattern(input_grid)

# Print the result in the required format
print('\n'.join(' '.join(str(cell) for cell in row) for row in output_grid))