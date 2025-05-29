def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = [row.split() for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    output = [row[:] for row in grid]
    
    # Process each row
    for r in range(rows):
        # If there's a 3 at the right edge
        if grid[r][cols-1] == '3':
            # Find all 1s positions in the row
            one_positions = [-1] + [i for i, x in enumerate(grid[r]) if x == '1'] + [cols]
            
            # Process segments between 1s from right to left
            for i in range(len(one_positions)-1, 0, -1):
                start, end = one_positions[i-1], one_positions[i]
                # If this segment contains 7s, convert them to 3s
                if i % 2 == 1:  # Alternate segments
                    for c in range(start+1, end):
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