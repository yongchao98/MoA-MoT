def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output = [[7 for _ in range(cols)] for _ in range(rows)]
    
    # Find adjacent 9s
    for i in range(rows):
        for j in range(cols-1):
            if input_grid[i][j] == 9 and input_grid[i][j+1] == 9:
                # Top row of 9s
                for c in range(7):
                    if 0 <= j-1+c < cols and 0 <= i-1 < rows:
                        output[i-1][j-1+c] = 9
                
                # Three rows of paired 3s
                for r in range(3):
                    if 0 <= i+r < rows:
                        # Left pair of 3s
                        if 0 <= j-1 < cols:
                            output[i+r][j-1] = 3
                            output[i+r][j] = 3
                        # Right pair of 3s
                        if 0 <= j+3 < cols:
                            output[i+r][j+3] = 3
                            output[i+r][j+4] = 3
                
                # Bottom row of 9s
                for c in range(7):
                    if 0 <= j-1+c < cols and 0 <= i+3 < rows:
                        output[i+3][j-1+c] = 9

    return output