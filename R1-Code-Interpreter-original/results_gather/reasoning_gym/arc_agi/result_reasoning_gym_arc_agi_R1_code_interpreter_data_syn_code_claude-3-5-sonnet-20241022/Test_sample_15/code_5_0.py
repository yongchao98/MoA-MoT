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
    
    # Create reflections at fixed distances from center
    pattern_size = 2  # Distance from center to reflect
    
    # For each position around center
    for i in range(center_row - pattern_size, center_row + pattern_size + 1):
        for j in range(center_col - pattern_size, center_col + pattern_size + 1):
            if 0 <= i < 10 and 0 <= j < 10:
                # Calculate relative position to center
                rel_row = i - center_row
                rel_col = j - center_col
                
                # Only reflect at specific distances
                if abs(rel_row) == pattern_size or abs(rel_col) == pattern_size:
                    # Find the value to reflect (look towards center)
                    source_row = center_row + (rel_row//2 if abs(rel_row) == pattern_size else rel_row)
                    source_col = center_col + (rel_col//2 if abs(rel_col) == pattern_size else rel_col)
                    
                    if 0 <= source_row < 10 and 0 <= source_col < 10:
                        value = input_grid[source_row][source_col]
                        if value != 2:
                            output_grid[i][j] = value
    
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