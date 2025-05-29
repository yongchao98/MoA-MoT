def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = [row.split() for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    output = [row[:] for row in grid]
    
    # Process each row
    for r in range(rows):
        if grid[r][cols-1] == '3':
            # Find all 1s positions from right to left
            one_positions = []
            for c in range(cols-1, -1, -1):
                if grid[r][c] == '1':
                    one_positions.append(c)
            one_positions.append(-1)  # Add left boundary
            
            # Process segments between 1s
            convert = True  # Start with converting (from right edge)
            last_pos = cols
            
            for pos in one_positions:
                if convert:
                    # Convert all 7s to 3s in this segment
                    for c in range(pos + 1, last_pos):
                        if grid[r][c] == '7':
                            output[r][c] = '3'
                last_pos = pos
                convert = not convert  # Toggle between converting and not converting
    
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