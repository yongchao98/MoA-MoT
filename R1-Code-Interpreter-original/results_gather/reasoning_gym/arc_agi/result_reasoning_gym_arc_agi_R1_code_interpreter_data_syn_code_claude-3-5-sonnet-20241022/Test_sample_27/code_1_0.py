def find_9_blocks(grid):
    height = len(grid)
    width = len(grid[0])
    blocks = []
    
    # Find all 2x2 blocks of 9's from top to bottom, left to right
    for i in range(height-1):
        for j in range(width-1):
            if (grid[i][j] == 9 and grid[i][j+1] == 9 and
                grid[i+1][j] == 9 and grid[i+1][j+1] == 9):
                blocks.append((i, j))
    return blocks

def apply_pattern(grid):
    height = len(grid)
    width = len(grid[0])
    output = [[3 for _ in range(width)] for _ in range(height)]
    
    # First copy all 9's
    for i in range(height):
        for j in range(width):
            if grid[i][j] == 9:
                output[i][j] = 9
    
    # Process blocks in order
    blocks = find_9_blocks(grid)
    
    # For each block, apply the pattern immediately
    for block_i, block_j in blocks:
        # Define the pattern positions
        positions = [
            (block_i-1, block_j-1, 5),  # top-left
            (block_i-1, block_j+2, 4),  # top-right
            (block_i+2, block_j-1, 1),  # bottom-left
            (block_i+2, block_j+2, 8)   # bottom-right
        ]
        
        # Apply each position if valid
        for pos_i, pos_j, value in positions:
            if 0 <= pos_i < height and 0 <= pos_j < width:
                # Check if position is not part of a 9's block
                is_nine = False
                for bi, bj in blocks:
                    if bi <= pos_i <= bi+1 and bj <= pos_j <= bj+1:
                        is_nine = True
                        break
                
                if not is_nine:
                    output[pos_i][pos_j] = value
    
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