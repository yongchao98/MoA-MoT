def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output = [[7 for _ in range(cols)] for _ in range(rows)]
    
    # Find 2x2 pattern of 9s or adjacent 9s
    for i in range(rows-1):
        for j in range(cols-1):
            if input_grid[i][j] == 9 and input_grid[i][j+1] == 9:
                # Create the pattern:
                # Top row of 9s
                for c in range(7):
                    if j-2+c < cols:
                        output[i-2][j-2+c] = 9
                
                # Middle rows with 3s
                for r in range(3):
                    if i-1+r < rows:
                        # Left pair of 3s
                        output[i-1+r][j-2] = 3
                        output[i-1+r][j-1] = 3
                        # Right pair of 3s
                        output[i-1+r][j+2] = 3
                        output[i-1+r][j+3] = 3
                
                # Bottom row of 9s
                for c in range(7):
                    if j-2+c < cols:
                        output[i+2][j-2+c] = 9

    return output

# Test input
test_input = [
    [7]*16 for _ in range(10)] + \
    [[7,7,7,7,7,7,7,9,9,7,7,7,7,7,7,7]] + \
    [[7]*16 for _ in range(5)]

# Generate output
output = transform_grid(test_input)
print("Output grid:")
for row in output:
    print(" ".join(map(str, row)))