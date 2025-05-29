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

def is_valid_position(grid, i, j):
    return 0 <= i < len(grid) and 0 <= j < len(grid[0])

def apply_pattern(grid):
    height = len(grid)
    width = len(grid[0])
    output = [[3 for _ in range(width)] for _ in range(height)]
    
    # Copy all 9's
    for i in range(height):
        for j in range(width):
            if grid[i][j] == 9:
                output[i][j] = 9
    
    blocks = find_9_blocks(grid)
    
    # Process each block
    for block_i, block_j in blocks:
        # Check each corner position
        corners = [
            (block_i-1, block_j-1, 5),  # top-left
            (block_i-1, block_j+2, 4),  # top-right
            (block_i+2, block_j-1, 1),  # bottom-left
            (block_i+2, block_j+2, 8)   # bottom-right
        ]
        
        for i, j, value in corners:
            if is_valid_position(grid, i, j):
                # Check if the position is not part of another 9's block
                is_valid = True
                for bi, bj in blocks:
                    if (bi <= i <= bi+1 and bj <= j <= bj+1):
                        is_valid = False
                        break
                if is_valid:
                    output[i][j] = value
    
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