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

def apply_pattern(grid):
    height = len(grid)
    width = len(grid[0])
    output = [[3 for _ in range(width)] for _ in range(height)]
    
    # First copy all 9's
    for i in range(height):
        for j in range(width):
            if grid[i][j] == 9:
                output[i][j] = 9
    
    # Process each 2x2 block of 9's
    blocks = find_9_blocks(grid)
    processed = set()  # Keep track of processed positions
    
    for block_i, block_j in blocks:
        # Define the four corner positions and their values
        corners = [
            (block_i-1, block_j-1, 5),  # top-left
            (block_i-1, block_j+2, 4),  # top-right
            (block_i+2, block_j-1, 1),  # bottom-left
            (block_i+2, block_j+2, 8)   # bottom-right
        ]
        
        for i, j, value in corners:
            # Check if position is valid and not already processed
            if (0 <= i < height and 0 <= j < width and 
                output[i][j] == 3 and  # position still has default value
                (i, j) not in processed):
                # Check if position is not part of any 9's block
                is_valid = True
                for bi, bj in blocks:
                    if (bi <= i <= bi+1 and bj <= j <= bj+1):
                        is_valid = False
                        break
                if is_valid:
                    output[i][j] = value
                    processed.add((i, j))
    
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