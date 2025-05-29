def generate_output(input_grid):
    n = len(input_grid)
    m = len(input_grid[0])
    output = [[1]*m for _ in range(n)]
    
    # Find special number in last row
    special_pos = next(i for i, x in enumerate(input_grid[-1]) if x != 1)
    special_num = input_grid[-1][special_pos]
    
    # Fill alternating columns from special position leftward
    for col in range(special_pos, -1, -2):
        for row in range(n):
            output[row][col] = special_num
    
    # Add 5s pattern
    # Start from special position and move leftward
    # First row: P-1, P-5, P-9, ...
    for pos in range(special_pos-1, -1, -4):
        output[0][pos] = 5
    
    # Last row: P-2, P-6, P-10, ...
    for pos in range(special_pos-2, -1, -4):
        output[-1][pos] = 5
    
    # Print result
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
test_input = [[1]*10 for _ in range(9)] + [[1,1,1,1,1,1,1,3,1,1]]
generate_output(test_input)