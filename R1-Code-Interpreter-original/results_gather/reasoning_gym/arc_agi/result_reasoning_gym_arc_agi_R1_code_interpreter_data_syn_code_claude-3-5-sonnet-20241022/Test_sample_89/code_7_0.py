def find_block_of_fives(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    def check_block(i, j):
        # Find maximum possible block of 5s starting at (i,j)
        if grid[i][j] != 5:
            return 0, 0
            
        # Find width
        width = 0
        for c in range(j, cols):
            if grid[i][c] != 5:
                break
            width += 1
            
        if width < 2:  # Need at least width of 2
            return 0, 0
            
        # Find height
        height = 0
        for r in range(i, rows):
            all_fives = True
            for c in range(j, j + width):
                if grid[r][c] != 5:
                    all_fives = False
                    break
            if not all_fives:
                break
            height += 1
            
        if height < 2:  # Need at least height of 2
            return 0, 0
            
        return height, width
    
    # Search for first valid block
    for i in range(rows-1):  # -1 because we need at least 2 rows
        for j in range(cols-1):  # -1 because we need at least 2 columns
            height, width = check_block(i, j)
            if height >= 2 and width >= 2:
                return (i, j, height, width)
                
    return None

def solve(input_grid):
    # Convert input string to 2D list
    grid = [[int(num) for num in row.split()] for row in input_grid.strip().split('\n')]
    
    # Create output grid as a copy of input grid
    output_grid = [row[:] for row in grid]
    
    # Find the first valid block
    block = find_block_of_fives(grid)
    
    # Replace the block with 7s if found
    if block:
        start_row, start_col, height, width = block
        for r in range(start_row, start_row + height):
            for c in range(start_col, start_col + width):
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