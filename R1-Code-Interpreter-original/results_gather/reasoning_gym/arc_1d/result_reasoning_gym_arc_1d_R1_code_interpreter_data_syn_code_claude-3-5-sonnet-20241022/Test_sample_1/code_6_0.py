def transform_grid(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Perform circular shift right by 2 positions
    n = len(grid)
    shift = 2
    output = grid[-shift:] + grid[:-shift]
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 6 4 4 7 8 8 6 9 8 5 8 2 1 8 0 0 0 0 0 0 0 0"
print(transform_grid(test_input))