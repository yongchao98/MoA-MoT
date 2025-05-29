def find_2x2_blocks(grid):
    blocks = []
    for i in range(len(grid)-1):
        for j in range(len(grid[0])-1):
            if grid[i][j] == 2 and grid[i][j+1] == 2 and grid[i+1][j] == 2 and grid[i+1][j+1] == 2:
                blocks.append((i, j))
    return blocks

def find_8s(grid):
    eights = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == 8:
                eights.append((i, j))
    return eights

def process_grid(input_grid):
    # Convert input string to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    output = [[0 for _ in range(10)] for _ in range(10)]
    
    # Find and copy 2x2 blocks to output
    blocks = find_2x2_blocks(grid)
    for block_i, block_j in blocks:
        for i in range(2):
            for j in range(2):
                output[block_i+i][block_j+j] = 2
    
    # Find original 8s
    eights = find_8s(grid)
    num_8s = len(eights)
    
    # For each block, place 8s in specific positions
    for block_i, block_j in blocks:
        positions = [
            (block_i, block_j-1),    # left side
            (block_i-1, block_j),    # above left
            (block_i-1, block_j+1),  # above right
            (block_i, block_j+2),    # right side
        ]
        
        # Place 8s in the defined positions
        for i in range(min(num_8s, 4)):
            pos = positions[i]
            if 0 <= pos[0] < 10 and 0 <= pos[1] < 10:
                output[pos[0]][pos[1]] = 8
    
    # Print the output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = """0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 8 0
0 0 0 8 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 2 2 0 0 0 0 0 0
8 0 2 2 0 0 0 0 0 8
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0"""

process_grid(test_input)