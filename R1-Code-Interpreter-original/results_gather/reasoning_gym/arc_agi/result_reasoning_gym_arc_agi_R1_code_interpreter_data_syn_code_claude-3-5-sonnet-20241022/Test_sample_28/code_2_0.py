def find_output(input_grid):
    # Find all non-6 numbers and their positions
    numbers = []
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 6:
                numbers.append((i, j, input_grid[i][j]))
    
    # Initialize output grid
    output = [[6 for _ in range(3)] for _ in range(3)]
    
    # Map based on position in input grid
    for i, j, val in numbers:
        # Top third of input (rows 0-3)
        if i < 3:
            if j < 3:
                output[0][0] = val
            elif j < 7:
                output[0][1] = val
            else:
                output[0][2] = val
        # Middle third of input (rows 3-6)
        elif i < 7:
            if j < 3:
                output[1][0] = val
            elif j < 7:
                output[1][1] = val
            else:
                output[1][2] = val
        # Bottom third of input (rows 7-9)
        else:
            if j < 3:
                output[2][0] = val
            elif j < 7:
                output[2][1] = val
            else:
                output[2][2] = val
    
    # Print output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [6, 6, 0, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 4, 6, 6, 6],
    [3, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 1, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 0, 6],
    [6, 5, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 3, 6],
    [6, 6, 6, 6, 6, 2, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 9, 6, 6]
]

find_output(test_input)