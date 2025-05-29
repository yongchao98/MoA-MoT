def transform_grid(input_grid):
    output_grid = [[2]*10 for _ in range(10)]
    
    # First, copy the input grid
    for i in range(10):
        for j in range(10):
            output_grid[i][j] = input_grid[i][j]
    
    # Find pattern boundaries
    min_row = max_row = min_col = max_col = None
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                if min_row is None:
                    min_row = max_row = i
                    min_col = max_col = j
                else:
                    min_row = min(min_row, i)
                    max_row = max(max_row, i)
                    min_col = min(min_col, j)
                    max_col = max(max_col, j)
    
    if min_row is None:
        return output_grid
        
    # Calculate center
    center_row = (min_row + max_row) // 2
    center_col = (min_col + max_col) // 2
    
    # For each non-2 value, create reflections
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                value = input_grid[i][j]
                
                # Calculate relative position to center
                rel_row = i - center_row
                rel_col = j - center_col
                
                # Create reflections
                positions = [
                    (center_row + rel_row, center_col + rel_col),    # original
                    (center_row + rel_row, center_col - rel_col),    # horizontal
                    (center_row - rel_row, center_col + rel_col),    # vertical
                    (center_row - rel_row, center_col - rel_col)     # diagonal
                ]
                
                # Apply all reflections
                for pos_row, pos_col in positions:
                    if 0 <= pos_row < 10 and 0 <= pos_col < 10:
                        # Don't overwrite original pattern
                        if input_grid[pos_row][pos_col] == 2:
                            output_grid[pos_row][pos_col] = value
    
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