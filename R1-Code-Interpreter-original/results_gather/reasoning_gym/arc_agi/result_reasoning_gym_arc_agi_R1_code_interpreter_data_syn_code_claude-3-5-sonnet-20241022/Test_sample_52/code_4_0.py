def process_grid(input_grid):
    n = len(input_grid)
    output = [[8 for _ in range(n)] for _ in range(n)]
    
    # Copy 2s from input to output
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 2:
                output[i][j] = 2
    
    # Process 4s first - create cross patterns
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 4:
                # Vertical line
                for r in range(max(0, i-3), min(n, i+4)):
                    if output[r][j] == 8:
                        output[r][j] = 4
                # Horizontal line
                for c in range(max(0, j-3), min(n, j+4)):
                    if output[i][c] == 8:
                        output[i][c] = 4
    
    # Connect special points in same row or column
    for i1 in range(n):
        for j1 in range(n):
            if input_grid[i1][j1] in [2, 4]:
                # Check same row
                for j2 in range(max(0, j1-4), min(n, j1+5)):
                    if j2 != j1 and input_grid[i1][j2] in [2, 4]:
                        for j in range(min(j1, j2), max(j1, j2)+1):
                            if output[i1][j] == 8:
                                output[i1][j] = 4
                
                # Check same column
                for i2 in range(max(0, i1-4), min(n, i1+5)):
                    if i2 != i1 and input_grid[i2][j1] in [2, 4]:
                        for i in range(min(i1, i2), max(i1, i2)+1):
                            if output[i][j1] == 8:
                                output[i][j1] = 4

    # Print the result
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
test_input = [
    [8,8,8,8,8,8,2,8,8,8,8],
    [8,2,8,8,8,8,8,8,8,8,2],
    [2,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,2,8,8,4,8,2,8,8],
    [8,8,8,8,8,8,8,8,8,8,8],
    [8,2,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8],
    [8,2,8,8,8,8,2,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,2]
]

process_grid(test_input)