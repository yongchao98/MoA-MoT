def find_9_blocks(grid):
    height = len(grid)
    width = len(grid[0])
    blocks = []
    
    for i in range(height-1):
        for j in range(width-1):
            if (grid[i][j] == 9 and grid[i][j+1] == 9 and
                grid[i+1][j] == 9 and grid[i+1][j+1] == 9):
                blocks.append((i, j))
    return sorted(blocks) # Sort blocks to ensure consistent processing

def apply_pattern(grid):
    height = len(grid)
    width = len(grid[0])
    output = [[3 for _ in range(width)] for _ in range(height)]
    
    # Copy all 9's first
    for i in range(height):
        for j in range(width):
            if grid[i][j] == 9:
                output[i][j] = 9
    
    # Find all 2x2 blocks of 9's
    blocks = find_9_blocks(grid)
    
    # For each block, apply the pattern
    for block_i, block_j in blocks:
        # Check if we can place all four numbers
        if (block_i > 0 and block_j > 0 and 
            block_i + 2 < height and block_j + 2 < width):
            
            # Check if positions are not part of another 9's block
            positions = [
                (block_i-1, block_j-1, 5),  # top-left
                (block_i-1, block_j+2, 4),  # top-right
                (block_i+2, block_j-1, 1),  # bottom-left
                (block_i+2, block_j+2, 8)   # bottom-right
            ]
            
            for pos_i, pos_j, value in positions:
                # Only place number if position is not part of a 9's block
                is_valid = True
                for bi, bj in blocks:
                    if (bi <= pos_i <= bi+1 and bj <= pos_j <= bj+1):
                        is_valid = False
                        break
                if is_valid and output[pos_i][pos_j] == 3:
                    output[pos_i][pos_j] = value
    
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

# Print the result
print('\n'.join(' '.join(str(cell) for cell in row) for row in output_grid))