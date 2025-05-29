def find_9_blocks(grid):
    height = len(grid)
    width = len(grid[0])
    blocks = []
    
    for i in range(height-1):
        for j in range(width-1):
            if (grid[i][j] == 9 and grid[i][j+1] == 9 and
                grid[i+1][j] == 9 and grid[i+1][j+1] == 9):
                blocks.append((i, j))
    return blocks

def is_valid_position(grid, i, j, blocks):
    height = len(grid)
    width = len(grid[0])
    
    # Check boundaries
    if not (0 <= i < height and 0 <= j < width):
        return False
        
    # Check if position is part of any 9's block
    for bi, bj in blocks:
        if (bi <= i <= bi+1 and bj <= j <= bj+1):
            return False
            
    return True

def apply_pattern(grid):
    height = len(grid)
    width = len(grid[0])
    output = [[3 for _ in range(width)] for _ in range(height)]
    
    # Copy 9's
    for i in range(height):
        for j in range(width):
            if grid[i][j] == 9:
                output[i][j] = 9
    
    # Find blocks
    blocks = find_9_blocks(grid)
    
    # Process each block
    for block_i, block_j in blocks:
        # Define pattern positions
        pattern = [
            (block_i-1, block_j-1, 5),  # top-left
            (block_i-1, block_j+2, 4),  # top-right
            (block_i+2, block_j-1, 1),  # bottom-left
            (block_i+2, block_j+2, 8)   # bottom-right
        ]
        
        # Apply pattern
        for i, j, value in pattern:
            if is_valid_position(grid, i, j, blocks):
                output[i][j] = value
    
    return output

# Test input
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

# Get output
output_grid = apply_pattern(input_grid)

# Print result
print('\n'.join(' '.join(str(cell) for cell in row) for row in output_grid))