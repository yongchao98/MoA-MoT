def find_block(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    # Check each possible starting position
    for i in range(rows-1):  # -1 because we need at least 2 rows
        for j in range(cols-1):  # -1 because we need at least 2 columns
            # Check if we have a valid block starting here
            if grid[i][j] not in [4, 5]:
                continue
                
            # Try to find the largest rectangular block of the same number
            val = grid[i][j]
            max_width = 1
            max_height = 1
            
            # Find width
            for w in range(j+1, cols):
                if grid[i][w] != val:
                    break
                max_width += 1
                
            # Find height
            valid_block = True
            for h in range(i+1, rows):
                for w in range(j, j + max_width):
                    if grid[h][w] != val:
                        valid_block = False
                        break
                if not valid_block:
                    break
                max_height += 1
            
            # Check if block is at least 2x2
            if max_height >= 2 and max_width >= 2:
                return (i, j, max_height, max_width)
    
    return None

def solve(input_grid):
    # Convert input string to 2D list
    grid = [[int(num) for num in row.split()] for row in input_grid.strip().split('\n')]
    
    # Create output grid as a copy of input grid
    output_grid = [row[:] for row in grid]
    
    # Find the first valid block
    block = find_block(grid)
    
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