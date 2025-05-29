def find_distinct_pattern(grid):
    # Find the distinct pattern in the grid
    distinct_values = set()
    for row in grid:
        for value in row:
            if value != 2:  # Assuming 2 is the background value
                distinct_values.add(value)
    
    # Extract the subgrid with the distinct pattern
    min_row, max_row = len(grid), 0
    min_col, max_col = len(grid[0]), 0
    
    for r in range(len(grid)):
        for c in range(len(grid[0])):
            if grid[r][c] in distinct_values:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)
    
    # Create the output grid
    output_grid = []
    for r in range(min_row, max_row + 1):
        output_grid.append(grid[r][min_col:max_col + 1])
    
    return output_grid

# Test input grid
test_input = [
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 8, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 8, 8, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 8, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [8, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
]

output = find_distinct_pattern(test_input)
for row in output:
    print(' '.join(map(str, row)))