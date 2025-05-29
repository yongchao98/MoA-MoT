def transform_grid(input_str):
    # Convert input string to list of integers
    grid = [int(x) for x in input_str.split()]
    
    # Create output grid
    output = [0] * len(grid)
    
    # For each position except the last, copy the next number
    for i in range(len(grid)-1):
        output[i] = grid[i+1]
    
    print(' '.join(map(str, output)))

test_input = "0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0"
transform_grid(test_input)