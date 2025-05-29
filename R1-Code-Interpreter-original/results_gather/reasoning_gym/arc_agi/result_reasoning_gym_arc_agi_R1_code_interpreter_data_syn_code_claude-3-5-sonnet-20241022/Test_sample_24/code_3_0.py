def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = [row.split() for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    output = [row[:] for row in grid]
    
    # First pass: horizontal spread from right edge
    three_rows = []
    for r in range(rows):
        if grid[r][cols-1] == '3':
            three_rows.append(r)
            # Find rightmost 1
            rightmost_1 = -1
            for c in range(cols-1, -1, -1):
                if grid[r][c] == '1':
                    rightmost_1 = c
                    break
            # Convert 7s to 3s from rightmost 1 to right edge
            for c in range(rightmost_1 + 1, cols):
                if grid[r][c] == '7':
                    output[r][c] = '3'
    
    # Second pass: vertical spread between rows with 3s
    if len(three_rows) > 1:
        for i in range(len(three_rows)-1):
            start_row = three_rows[i]
            end_row = three_rows[i+1]
            # Check if there's a blocking row of 1s between
            blocked = False
            for r in range(start_row+1, end_row):
                if '1' in grid[r]:
                    blocked = True
                    break
            if not blocked:
                # Fill all 7s between these rows with 3s
                for r in range(start_row+1, end_row):
                    for c in range(cols):
                        if grid[r][c] == '7':
                            output[r][c] = '3'
    
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