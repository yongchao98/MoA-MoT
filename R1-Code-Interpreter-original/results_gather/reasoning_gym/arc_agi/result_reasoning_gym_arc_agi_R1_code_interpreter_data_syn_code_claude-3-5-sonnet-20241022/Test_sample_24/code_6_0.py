def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = [row.split() for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    output = [row[:] for row in grid]
    
    # First pass: handle rows ending with 3
    three_rows = []
    for r in range(rows):
        if grid[r][cols-1] == '3':
            three_rows.append(r)
            # Find all 1s positions in the row
            one_positions = [-1] + [i for i, x in enumerate(grid[r]) if x == '1']
            one_positions.append(cols)
            
            # Process segments between 1s from right to left
            fill = True  # Start with filling (from right edge)
            for i in range(len(one_positions)-1, 0, -1):
                if fill:
                    # Fill all 7s with 3s in this segment
                    for c in range(one_positions[i-1]+1, one_positions[i]):
                        if grid[r][c] == '7':
                            output[r][c] = '3'
                fill = not fill  # Alternate between filling and not filling
    
    # Second pass: fill between rows ending with 3
    for i in range(len(three_rows)-1):
        start_row = three_rows[i]
        end_row = three_rows[i+1]
        # Fill all rows between if there are no 1s blocking
        can_fill = True
        for r in range(start_row+1, end_row):
            if '1' in grid[r]:
                can_fill = False
                break
        
        if can_fill:
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