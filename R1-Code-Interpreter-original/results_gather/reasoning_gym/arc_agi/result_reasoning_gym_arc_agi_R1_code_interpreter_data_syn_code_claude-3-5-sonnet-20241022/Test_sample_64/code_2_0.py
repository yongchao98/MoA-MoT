def transpose_grid(grid):
    # Create transposed grid
    transposed = [[grid[j][i] for j in range(3)] for i in range(3)]
    
    print("Input grid:")
    for row in grid:
        print(row)
    print("\nTransposed grid:")
    for row in transposed:
        print(row)
    return transposed

# Test with the given input
test_input = [[4,0,0], [0,4,4], [4,4,4]]
result = transpose_grid(test_input)