def transform_grid(input_grid):
    output_grid = [[2]*10 for _ in range(10)]
    
    # First, copy the input grid
    for i in range(10):
        for j in range(10):
            output_grid[i][j] = input_grid[i][j]
    
    # Find the pattern center
    center_found = False
    center_row = center_col = 0
    
    # Scan for non-2 values to find center
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                min_row = i
                center_found = True
                break
        if center_found:
            break
    
    # Find the extent of the pattern
    max_row = min_row
    min_col = max_col = 0
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                max_row = max(max_row, i)
                if min_col == 0 or j < min_col:
                    min_col = j
                max_col = max(max_col, j)
    
    # Calculate center
    center_row = (min_row + max_row) // 2
    center_col = (min_col + max_col) // 2
    
    # Apply the transformation
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                # Calculate distances from center
                row_dist = i - center_row
                col_dist = j - center_col
                
                # Mirror positions
                mirror_positions = [
                    (center_row + row_dist, center_col - col_dist),  # horizontal mirror
                    (center_row - row_dist, center_col + col_dist),  # vertical mirror
                    (center_row - row_dist, center_col - col_dist)   # diagonal mirror
                ]
                
                # Apply mirrors
                for mirror_row, mirror_col in mirror_positions:
                    if 0 <= mirror_row < 10 and 0 <= mirror_col < 10:
                        output_grid[mirror_row][mirror_col] = input_grid[i][j]
    
    return output_grid

# Test input
input_grid = [
    [2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,1,2,2,2,2],
    [2,2,2,2,0,1,0,0,0,2],
    [2,2,2,2,1,0,1,2,2,2],
    [2,2,2,2,0,1,0,0,0,2],
    [2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2]
]

result = transform_grid(input_grid)
for row in result:
    print(' '.join(map(str, row)))