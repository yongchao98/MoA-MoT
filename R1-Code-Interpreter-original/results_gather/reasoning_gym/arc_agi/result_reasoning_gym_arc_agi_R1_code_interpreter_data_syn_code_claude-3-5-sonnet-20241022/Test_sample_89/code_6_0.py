def find_valid_block(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    def check_valid_block(i, j, h, w):
        # Check if block contains only 4s and 5s
        for r in range(i, min(i + h, rows)):
            for c in range(j, min(j + w, cols)):
                if grid[r][c] not in [4, 5]:
                    return False
        return True
    
    # Try each possible starting position
    for i in range(rows-1):
        for j in range(cols-1):
            # Check 2x2 block first
            if check_valid_block(i, j, 2, 2):
                # Try to expand the block
                max_height = 2
                max_width = 2
                
                # Try expanding width
                while j + max_width < cols and check_valid_block(i, j, 2, max_width + 1):
                    max_width += 1
                
                # Try expanding height
                while i + max_height < rows and check_valid_block(i, j, max_height + 1, max_width):
                    max_height += 1
                
                return (i, j, max_height, max_width)
    return None

def solve(input_grid):
    # Convert input string to 2D list
    grid = [[int(num) for num in row.split()] for row in input_grid.strip().split('\n')]
    
    # Create output grid as a copy of input grid
    output_grid = [row[:] for row in grid]
    
    # Find the first valid block
    block = find_valid_block(grid)
    
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