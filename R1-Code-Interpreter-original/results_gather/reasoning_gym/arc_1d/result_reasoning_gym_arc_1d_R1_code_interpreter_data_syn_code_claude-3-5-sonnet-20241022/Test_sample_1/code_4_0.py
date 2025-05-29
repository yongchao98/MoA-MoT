def transform_grid(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Create output array
    output = [0] * 30
    
    # Place last two numbers at start
    output[0] = grid[-2]
    output[1] = grid[-1]
    
    # Place third and fourth last numbers at end
    output[-2] = grid[-4]
    output[-1] = grid[-3]
    
    # Shift everything else right by 2
    for i in range(len(grid)-4):
        output[i+2] = grid[i]
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 6 4 4 7 8 8 6 9 8 5 8 2 1 8 0 0 0 0 0 0 0 0"
print(transform_grid(test_input))