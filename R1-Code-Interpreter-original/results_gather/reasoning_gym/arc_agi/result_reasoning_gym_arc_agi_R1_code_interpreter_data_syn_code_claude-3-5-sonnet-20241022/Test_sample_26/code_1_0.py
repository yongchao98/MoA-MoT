def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    def find_pairs():
        pairs = []
        for i in range(rows-1):
            for j in range(cols-1):
                if input_grid[i][j] == input_grid[i][j+1] == input_grid[i+1][j] == input_grid[i+1][j+1] == 1:
                    pairs.append((i, j))
        return pairs
    
    pairs = find_pairs()
    
    # Check each pair
    for i, j in pairs:
        # Check if there's any pair above this one
        has_pair_above = False
        for y in range(0, i):
            if (y, j) in pairs:
                has_pair_above = True
                break
        
        # If there's a pair above, transform this pair to 0's
        if has_pair_above:
            output_grid[i][j] = output_grid[i][j+1] = output_grid[i+1][j] = output_grid[i+1][j+1] = 0
            
    # Print the result in one line
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