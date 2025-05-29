def transform_grid(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    n = len(grid)
    
    # Initialize output array
    output = [0] * n
    
    # Step 1: Last two numbers go to front
    output[0] = grid[-2]
    output[1] = grid[-1]
    
    # Step 2: Shift remaining numbers right by 2
    for i in range(n-2):
        output[i+2] = grid[i]
    
    # Step 3: Last two positions are zeros
    output[-2] = 0
    output[-1] = 0
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 6 4 4 7 8 8 6 9 8 5 8 2 1 8 0 0 0 0 0 0 0 0"
print(transform_grid(test_input))