def transform_grid(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    n = len(grid)
    
    # Initialize output array
    output = [0] * n
    
    # Step 1: Move last two numbers to front
    output[0] = grid[n-2]
    output[1] = grid[n-1]
    
    # Step 2: Shift remaining numbers right by two
    for i in range(n-4):
        output[i+2] = grid[i]
    
    # Step 3: Move (n-4)th and (n-3)th numbers to end
    output[n-2] = grid[n-4]
    output[n-1] = grid[n-3]
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 6 4 4 7 8 8 6 9 8 5 8 2 1 8 0 0 0 0 0 0 0 0"
print(transform_grid(test_input))