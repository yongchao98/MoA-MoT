def shift_sequence(input_grid):
    # Convert to string for easier manipulation
    grid = [str(x) for x in input_grid]
    
    # Find the first and last non-zero number
    start = 0
    end = 0
    for i in range(len(grid)):
        if grid[i] != '0':
            start = i
            break
    for i in range(len(grid)-1, -1, -1):
        if grid[i] != '0':
            end = i
            break
    
    # Get the number being used
    num = grid[start]
    length = end - start + 1
    
    # Create output grid
    output = ['0'] * 25
    new_start = start + 2  # Shift right by 2
    for i in range(length):
        output[new_start + i] = num
    
    # Print the result
    print(' '.join(output))

# Test input
test = [0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0]
shift_sequence(test)