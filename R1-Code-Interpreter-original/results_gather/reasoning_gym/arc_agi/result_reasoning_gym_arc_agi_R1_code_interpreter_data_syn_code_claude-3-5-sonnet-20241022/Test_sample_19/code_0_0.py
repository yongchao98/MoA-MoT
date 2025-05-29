def find_2x2_blocks(grid):
    blocks = []
    for i in range(len(grid)-1):
        for j in range(len(grid[0])-1):
            if grid[i][j] == 2 and grid[i][j+1] == 2 and grid[i+1][j] == 2 and grid[i+1][j+1] == 2:
                blocks.append((i, j))  # Store top-left corner
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
    
    # Find 2x2 blocks and copy them to output
    blocks = find_2x2_blocks(grid)
    for block_i, block_j in blocks:
        for i in range(2):
            for j in range(2):
                output[block_i+i][block_j+j] = 2
    
    # Find 8s and process them
    eights = find_8s(grid)
    for eight_i, eight_j in eights:
        # For each 8, find closest 2x2 block
        if not blocks:
            continue
            
        min_dist = float('inf')
        new_pos = None
        
        for block_i, block_j in blocks:
            # Try positions at same level and one row above block
            possible_positions = [
                (block_i, block_j-1), (block_i, block_j+2),  # sides
                (block_i-1, block_j-1), (block_i-1, block_j), 
                (block_i-1, block_j+1), (block_i-1, block_j+2)  # above
            ]
            
            for pos_i, pos_j in possible_positions:
                if 0 <= pos_i < 10 and 0 <= pos_j < 10:
                    dist = abs(eight_i - block_i) + abs(eight_j - block_j)
                    if dist < min_dist:
                        min_dist = dist
                        new_pos = (pos_i, pos_j)
        
        if new_pos:
            output[new_pos[0]][new_pos[1]] = 8

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