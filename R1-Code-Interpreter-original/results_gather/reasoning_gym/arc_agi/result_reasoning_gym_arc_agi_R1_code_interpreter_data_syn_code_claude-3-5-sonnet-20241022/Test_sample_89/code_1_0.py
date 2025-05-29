def find_largest_block(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    def check_block(i, j, target):
        # Find maximum block of target number (4 or 5)
        if grid[i][j] != target:
            return 0, 0
            
        # Find maximum width
        width = 0
        for c in range(j, cols):
            if grid[i][c] != target:
                break
            width += 1
            
        if width < 2:
            return 0, 0
            
        # Find maximum height
        height = 0
        for r in range(i, rows):
            all_same = True
            for c in range(j, j + width):
                if grid[r][c] != target:
                    all_same = False
                    break
            if not all_same:
                break
            height += 1
            
        if height < 2:
            return 0, 0
            
        return height, width
    
    # Find largest valid block
    max_size = 0
    best_block = None
    
    for i in range(rows-1):
        for j in range(cols-1):
            # Try both 4s and 5s
            for target in [4, 5]:
                height, width = check_block(i, j, target)
                if height >= 2 and width >= 2:
                    size = height * width
                    if size > max_size:
                        max_size = size
                        best_block = (i, j, height, width)
                        
    return best_block

def solve(input_grid):
    # Convert input string to 2D list
    grid = [[int(num) for num in row.split()] for row in input_grid.strip().split('\n')]
    
    # Create output grid as a copy of input grid
    output_grid = [row[:] for row in grid]
    
    # Find the largest valid block
    block = find_largest_block(grid)
    
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