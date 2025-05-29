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

def get_valid_positions(block_i, block_j):
    # Get all possible positions around the 2x2 block
    positions = []
    # Positions above the block
    for j in range(block_j-1, block_j+3):
        positions.append((block_i-1, j))
    # Positions at the same level
    positions.append((block_i, block_j-1))
    positions.append((block_i, block_j+2))
    # Positions below the block
    positions.append((block_i+2, block_j-1))
    positions.append((block_i+2, block_j+2))
    return [(i, j) for i, j in positions if 0 <= i < 10 and 0 <= j < 10]

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
    
    # Process 8s
    eights = find_8s(grid)
    used_positions = set()
    
    for block_i, block_j in blocks:
        valid_positions = get_valid_positions(block_i, block_j)
        local_eights = eights[:3]  # Take up to 3 8s per block
        eights = eights[3:]  # Remove used 8s
        
        # Distribute 8s around the block
        for _ in range(len(local_eights)):
            if not valid_positions:
                break
            # Choose position preferring those above the block
            for pos in valid_positions:
                if pos not in used_positions:
                    output[pos[0]][pos[1]] = 8
                    used_positions.add(pos)
                    valid_positions.remove(pos)
                    break
    
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