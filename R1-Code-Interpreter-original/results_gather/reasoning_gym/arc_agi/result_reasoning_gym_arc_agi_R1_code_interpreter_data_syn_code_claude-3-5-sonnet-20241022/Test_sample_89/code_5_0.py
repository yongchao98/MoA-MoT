def find_consecutive_block(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    def check_block(i, j, val):
        # Try to find maximum block size of consecutive same numbers
        max_width = 1
        while j + max_width < cols and grid[i][j + max_width] == val:
            max_width += 1
            
        max_height = 1
        for h in range(i + 1, rows):
            if all(grid[h][c] == val for c in range(j, j + max_width)):
                max_height += 1
            else:
                break
                
        return max_height, max_width
    
    # Search for first valid block
    for i in range(rows-1):
        for j in range(cols-1):
            val = grid[i][j]
            if val not in [4, 5]:
                continue
                
            height, width = check_block(i, j, val)
            if height >= 2 and width >= 2:
                # Verify it's the largest possible block
                return (i, j, height, width)
    
    return None

def solve(input_grid):
    # Convert input string to 2D list
    grid = [[int(num) for num in row.split()] for row in input_grid.strip().split('\n')]
    
    # Create output grid as a copy of input grid
    output_grid = [row[:] for row in grid]
    
    # Find the first consecutive block
    block = find_consecutive_block(grid)
    
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