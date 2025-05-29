def find_identical_block(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    def check_block(i, j, val):
        # Try to find maximum block size starting at (i,j)
        max_width = 1
        max_height = 1
        
        # Find maximum width
        for w in range(j + 1, cols):
            if grid[i][w] != val:
                break
            max_width += 1
            
        # Find maximum height
        for h in range(i + 1, rows):
            all_match = True
            for w in range(j, j + max_width):
                if grid[h][w] != val:
                    all_match = False
                    break
            if not all_match:
                break
            max_height += 1
            
        return max_height, max_width if max_height >= 2 and max_width >= 2 else (0, 0)
    
    # Search for first valid block
    for i in range(rows):
        for j in range(cols):
            val = grid[i][j]
            if val in [4, 5]:
                height, width = check_block(i, j, val)
                if height >= 2 and width >= 2:
                    return (i, j, height, width)
    
    return None

def solve(input_grid):
    # Convert input string to 2D list
    grid = [[int(num) for num in row.split()] for row in input_grid.strip().split('\n')]
    
    # Create output grid as a copy of input grid
    output_grid = [row[:] for row in grid]
    
    # Find the first identical block
    block = find_identical_block(grid)
    
    # Replace the block with 7s if found
    if block:
        i, j, height, width = block
        for r in range(i, i + height):
            for c in range(j, j + width):
                output_grid[r][c] = 7
    
    # Convert back to string format
    return '\n'.join(' '.join(str(num) for num in row) for row in output_grid)

# Test input
test_input = """4 4 5 5
5 4 4 4
5 4 4 5
5 4 4 4
4 5 4 4
4 5 4 5
5 4 5 5
4 4 4 5
4 5 4 5
5 5 5 4
4 5 4 5
4 4 4 5
4 5 5 5
4 4 5 5
4 4 5 5
5 5 5 4
5 4 4 5
5 4 4 4
5 4 4 5
5 5 4 4
4 5 5 5
4 5 5 5
5 5 5 5
5 5 4 5"""

print(solve(test_input))