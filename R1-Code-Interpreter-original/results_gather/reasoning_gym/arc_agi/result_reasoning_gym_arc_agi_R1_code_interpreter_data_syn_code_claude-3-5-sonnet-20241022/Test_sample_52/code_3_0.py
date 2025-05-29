def process_grid(input_grid):
    n = len(input_grid)
    output = [[8 for _ in range(n)] for _ in range(n)]
    
    # First, copy all 2s from input to output
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 2:
                output[i][j] = 2
    
    # Process 4s first
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] == 4:
                # Create cross pattern
                for r in range(max(0, i-3), min(n, i+4)):
                    if output[r][j] == 8:
                        output[r][j] = 4
                for c in range(max(0, j-3), min(n, j+4)):
                    if output[i][c] == 8:
                        output[i][c] = 4
    
    # Connect nearby 2s and 4s
    special_points = [(i, j) for i in range(n) for j in range(n) 
                     if input_grid[i][j] in [2, 4]]
    
    for i1, j1 in special_points:
        for i2, j2 in special_points:
            if (i1, j1) != (i2, j2):
                # If points are within 4 steps of each other
                if abs(i1 - i2) <= 4 and abs(j1 - j2) <= 4:
                    # Create rectangular path
                    if i1 == i2 or j1 == j2:  # If in same row or column
                        if i1 == i2:  # Same row
                            for j in range(min(j1, j2), max(j1, j2) + 1):
                                if output[i1][j] == 8:
                                    output[i1][j] = 4
                        else:  # Same column
                            for i in range(min(i1, i2), max(i1, i2) + 1):
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