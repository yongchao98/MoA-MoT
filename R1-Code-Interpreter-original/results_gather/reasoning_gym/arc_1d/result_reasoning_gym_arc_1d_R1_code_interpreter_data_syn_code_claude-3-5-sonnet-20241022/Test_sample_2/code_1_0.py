def solve_grid(input_grid):
    # Convert input string to list of integers
    grid = [int(x) for x in input_grid.split()]
    output = grid.copy()
    
    # Find first non-zero number
    first_non_zero = None
    first_non_zero_pos = None
    for i in range(len(grid)):
        if grid[i] != 0:
            first_non_zero = grid[i]
            first_non_zero_pos = i
            break
    
    # If found, fill all leading zeros with this number
    if first_non_zero is not None:
        for i in range(first_non_zero_pos):
            output[i] = first_non_zero
    
    # Print result
    print(' '.join(map(str, output)))

# Test with the given input
test_input = "0 0 3 3 3 3 3 3 3 3 3 3 3 3 2 0 0 0 0 0 0"
solve_grid(test_input)