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
    
    # For each special point, create paths to nearest special points
    for i, j in special_points:
        # Create a cross pattern around 4s
        if input_grid[i][j] == 4:
            # Vertical line
            for r in range(i-4, i+5):
                if 0 <= r < n:
                    output[r][j] = 4
            # Horizontal line
            for c in range(j-4, j+5):
                if 0 <= c < n:
                    output[i][c] = 4
        # Create paths between special points
        for i2, j2 in special_points:
            if (i, j) != (i2, j2):
                # If points are close enough, connect them
                if abs(i - i2) <= 4 and abs(j - j2) <= 4:
                    # Fill the rectangle between points with 4s
                    min_i, max_i = min(i, i2), max(i, i2)
                    min_j, max_j = min(j, j2), max(j, j2)
                    for r in range(min_i, max_i + 1):
                        for c in range(min_j, max_j + 1):
                            if output[r][c] == 8:  # Don't overwrite 2s
                                output[r][c] = 4

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