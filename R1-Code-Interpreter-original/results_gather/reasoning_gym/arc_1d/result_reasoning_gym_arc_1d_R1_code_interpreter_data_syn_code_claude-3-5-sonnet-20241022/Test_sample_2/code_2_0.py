def solve_grid(input_grid):
    # Convert input string to list of integers
    grid = [int(x) for x in input_grid.split()]
    output = grid.copy()
    
    # Find the last different number and its position
    last_num = None
    last_pos = None
    for i in range(len(grid)-1, -1, -1):
        if grid[i] != 0:
            if last_num is None or last_num == grid[i]:
                last_num = grid[i]
                last_pos = i
            else:
                break
    
    # Fill trailing zeros after last different number
    if last_num is not None:
        for i in range(last_pos + 1, len(grid)):
            if grid[i] == 0:
                output[i] = last_num
    
    # Print result
    print(' '.join(map(str, output)))

# Test with the given input
test_input = "0 0 3 3 3 3 3 3 3 3 3 3 3 3 2 0 0 0 0 0 0"
solve_grid(test_input)