def find_mixed_block(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    # Check each possible starting position
    for i in range(rows-1):
        for j in range(cols-1):
            # Try different block sizes starting from 2x2
            for height in range(2, min(rows-i+1, 5)):
                for width in range(2, min(cols-j+1, 5)):
                    # Check if this block contains both 4s and 5s
                    has_4 = False
                    has_5 = False
                    valid_block = True
                    
                    for r in range(i, min(i+height, rows)):
                        for c in range(j, min(j+width, cols)):
                            if grid[r][c] == 4:
                                has_4 = True
                            elif grid[r][c] == 5:
                                has_5 = True
                            else:
                                valid_block = False
                                break
                        if not valid_block:
                            break
                    
                    if valid_block and has_4 and has_5:
                        return (i, j, height, width)
    return None

def solve(input_grid):
    # Convert input string to 2D list
    grid = [[int(num) for num in row.split()] for row in input_grid.strip().split('\n')]
    
    # Create output grid as a copy of input grid
    output_grid = [row[:] for row in grid]
    
    # Find the first valid block
    block = find_mixed_block(grid)
    
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