def process_grid(input_grid):
    n = len(input_grid)
    output = [[8 for _ in range(n)] for _ in range(n)]
    
    # First, copy all 2s from input to output
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 2:
                output[i][j] = 2
    
    # Find all special points (2s and 4s)
    special_points = []
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] in [2, 4]:
                special_points.append((i, j))
    
    # Process each special point
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 4:
                # Create vertical path
                for r in range(i-3, i+4):
                    if 0 <= r < n and output[r][j] == 8:
                        output[r][j] = 4
                # Create horizontal path
                for c in range(j-3, j+4):
                    if 0 <= c < n and output[i][c] == 8:
                        output[i][c] = 4
            elif input_grid[i][j] == 2:
                # Find closest special point in same column or row
                for r in range(i-3, i+4):
                    if 0 <= r < n and r != i:
                        if input_grid[r][j] in [2, 4]:
                            # Fill vertical path
                            for k in range(min(r, i), max(r, i)+1):
                                if output[k][j] == 8:
                                    output[k][j] = 4
                for c in range(j-3, j+4):
                    if 0 <= c < n and c != j:
                        if input_grid[i][c] in [2, 4]:
                            # Fill horizontal path
                            for k in range(min(c, j), max(c, j)+1):
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