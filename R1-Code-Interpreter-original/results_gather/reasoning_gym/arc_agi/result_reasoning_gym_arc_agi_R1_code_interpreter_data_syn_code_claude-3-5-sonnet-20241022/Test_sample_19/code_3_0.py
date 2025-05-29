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

def manhattan_distance(p1, p2):
    return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])

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
    
    # Find original 8 positions
    original_8s = find_8s(grid)
    
    # For each block, define possible positions for 8s
    for block_i, block_j in blocks:
        positions = [
            (block_i, block_j-1),    # left side
            (block_i-1, block_j),    # above middle
            (block_i, block_j+2),    # right side
        ]
        
        # Sort original 8s by distance to this block
        sorted_8s = sorted(original_8s, 
                         key=lambda pos: manhattan_distance(pos, (block_i, block_j)))
        
        # Place 8s in the closest valid positions
        for eight in sorted_8s[:3]:  # Take up to 3 8s per block
            min_dist = float('inf')
            best_pos = None
            
            for pos in positions:
                if 0 <= pos[0] < 10 and 0 <= pos[1] < 10:
                    dist = manhattan_distance(eight, pos)
                    if dist < min_dist:
                        min_dist = dist
                        best_pos = pos
            
            if best_pos:
                output[best_pos[0]][best_pos[1]] = 8
                positions.remove(best_pos)
    
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