def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = [row.split() for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    output = [row[:] for row in grid]
    
    # Process each row that ends with 3, from bottom to top
    for r in range(rows-1, -1, -1):
        if grid[r][cols-1] == '3':
            # Handle the current row
            rightmost_1 = -1
            for c in range(cols-1, -1, -1):
                if grid[r][c] == '1':
                    rightmost_1 = c
                    break
            # Convert 7s to 3s from rightmost 1 to right edge
            for c in range(rightmost_1 + 1, cols):
                if grid[r][c] == '7':
                    output[r][c] = '3'
            
            # Spread upward until hitting a row with 1s
            for prev_row in range(r-1, -1, -1):
                if '1' in grid[prev_row]:
                    break
                for c in range(cols):
                    if grid[prev_row][c] == '7':
                        output[prev_row][c] = '3'
    
    return '\n'.join(' '.join(row) for row in output)

# Test input
test_input = """7 7 7 1 7 7 7 1 7 7 7 1 7
7 7 7 7 7 7 7 7 1 1 1 7 7
7 7 7 1 1 7 1 7 7 7 7 1 7
1 7 7 7 1 7 1 7 7 1 7 7 3
7 7 7 7 7 7 7 7 7 7 7 1 7
1 1 7 7 1 7 1 1 7 7 7 7 7
7 7 1 1 7 7 7 1 7 7 7 1 7
7 1 7 1 7 7 1 7 7 1 7 7 3
7 7 7 7 7 1 7 1 7 7 7 7 7
7 1 7 7 7 7 7 7 1 1 1 7 3
7 1 7 7 7 1 7 7 1 7 1 7 1
1 1 7 7 7 7 7 7 7 1 7 7 1
1 7 7 1 1 7 7 7 7 7 7 7 7"""

print(transform_grid(test_input))