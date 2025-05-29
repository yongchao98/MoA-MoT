def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    def is_pair(i, j):
        if i >= rows-1 or j >= cols-1:
            return False
        return (input_grid[i][j] == input_grid[i][j+1] == 
                input_grid[i+1][j] == input_grid[i+1][j+1] == 1)
    
    # Find all 2x2 blocks of 1's
    pairs = []
    for i in range(rows-1):
        for j in range(cols-1):
            if is_pair(i, j):
                pairs.append((i, j))
    
    # For each pair, check if there's any pair above it
    for i, j in pairs:
        for above_pair in pairs:
            if above_pair[1] == j and above_pair[0] < i:
                # Transform this pair to 0's
                output_grid[i][j] = output_grid[i][j+1] = 0
                output_grid[i+1][j] = output_grid[i+1][j+1] = 0
                break
    
    # Print the result
    for row in output_grid:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [8,8,8,8,8,1,1],
    [8,8,8,8,8,1,1],
    [8,8,1,1,8,8,8],
    [8,8,1,1,8,8,8],
    [8,8,8,8,1,1,8],
    [8,8,8,8,1,1,8],
    [8,8,8,8,8,8,8],
    [8,1,1,8,8,8,8],
    [8,1,1,8,8,1,1],
    [8,8,8,8,8,1,1],
    [8,8,8,8,8,8,8],
    [8,8,8,8,1,1,8],
    [8,8,8,8,1,1,8],
    [8,1,1,8,8,8,8],
    [8,1,1,8,8,8,8]
]

transform_grid(test_input)