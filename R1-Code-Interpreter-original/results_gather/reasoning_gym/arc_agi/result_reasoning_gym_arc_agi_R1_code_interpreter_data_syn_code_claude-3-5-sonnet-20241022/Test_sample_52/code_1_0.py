def process_grid(input_grid):
    n = len(input_grid)
    output = [[8 for _ in range(n)] for _ in range(n)]
    
    # First, copy all 2s from input to output
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 2:
                output[i][j] = 2
    
    # Process each 4 in the input
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 4:
                # Create diamond pattern around 4
                for di in range(-4, 5):
                    for dj in range(-4, 5):
                        if abs(di) + abs(dj) <= 4:  # Diamond shape
                            ni, nj = i + di, j + dj
                            if 0 <= ni < n and 0 <= nj < n:
                                if output[ni][nj] == 8:  # Don't overwrite 2s
                                    output[ni][nj] = 4
    
    # Process each 2 in the input
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 2:
                # Create vertical and horizontal lines of 4s
                for k in range(max(0, i-4), min(n, i+5)):
                    if output[k][j] == 8:
                        output[k][j] = 4
                for k in range(max(0, j-4), min(n, j+5)):
                    if output[i][k] == 8:
                        output[i][k] = 4

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