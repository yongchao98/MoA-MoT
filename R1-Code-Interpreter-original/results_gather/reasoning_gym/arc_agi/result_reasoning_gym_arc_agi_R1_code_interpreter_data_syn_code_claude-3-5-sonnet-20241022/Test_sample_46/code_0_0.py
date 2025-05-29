def find_blocks(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    blocks = []
    
    def get_block(r, c, val):
        if (r, c) in visited or r < 0 or r >= rows or c < 0 or c >= cols or grid[r][c] != val:
            return None
        
        # Find the extent of the block
        max_r, max_c = r, c
        for i in range(r, rows):
            if grid[i][c] != val:
                break
            max_r = i
        for j in range(c, cols):
            if grid[r][j] != val:
                break
            max_c = j
            
        # Verify if it's a valid rectangular block
        for i in range(r, max_r + 1):
            for j in range(c, max_c + 1):
                if grid[i][j] != val:
                    return None
                visited.add((i, j))
        
        return (r, c, max_r, max_c, val)
    
    # Find all blocks
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in visited and grid[i][j] != 7:
                block = get_block(i, j, grid[i][j])
                if block:
                    blocks.append(block)
    
    return blocks

def transform_grid(input_grid):
    rows, cols = len(input_grid), len(input_grid[0])
    output = [[7 for _ in range(cols)] for _ in range(rows)]
    
    blocks = find_blocks(input_grid)
    
    # Process each block
    for block in blocks:
        r1, c1, r2, c2, val = block
        # Check if block is large enough (at least 2x2 or 3x3)
        if (r2 - r1 + 1) * (c2 - c1 + 1) >= 4:
            for i in range(r1, r2 + 1):
                for j in range(c1, c2 + 1):
                    output[i][j] = val
    
    return output

# Test input
test_input = [
    [6,7,7,7,6,7,7,7,7,6],
    [7,7,6,6,6,7,7,7,7,7],
    [7,7,6,6,6,7,6,7,6,7],
    [7,7,6,6,6,7,7,7,7,7],
    [7,7,6,6,6,7,7,7,7,7],
    [6,7,6,6,6,7,7,7,7,7],
    [7,7,7,7,7,7,6,6,7,7],
    [7,6,7,7,7,7,6,6,7,7],
    [7,7,7,7,7,7,6,6,6,7],
    [7,7,6,7,7,6,7,7,7,7]
]

result = transform_grid(test_input)
for row in result:
    print(' '.join(map(str, row)))