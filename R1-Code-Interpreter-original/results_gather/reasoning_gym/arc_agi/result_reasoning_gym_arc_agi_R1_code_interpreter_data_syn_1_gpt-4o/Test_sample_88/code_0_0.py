def find_pairs(grid):
    rows_with_3 = set()
    cols_with_6 = set()
    
    # Find all rows with a 3
    for i, row in enumerate(grid):
        if 3 in row:
            rows_with_3.add(i)
    
    # Find all columns with a 6
    for j in range(len(grid[0])):
        for i in range(len(grid)):
            if grid[i][j] == 6:
                cols_with_6.add(j)
                break
    
    # Form pairs and sort them
    pairs = sorted((r, c) for r in rows_with_3 for c in cols_with_6)
    return pairs

# Test input grid
test_input = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 0, 6, 0, 0, 0, 3, 0, 0, 0],
    [0, 3, 0, 0, 0, 3, 0, 6, 0, 0, 0, 3, 0, 0, 0],
    [0, 3, 0, 0, 0, 3, 0, 6, 0, 0, 0, 3, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6, 0],
    [0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6, 0],
    [0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Find and print the output pairs
output_pairs = find_pairs(test_input)
print(output_pairs)