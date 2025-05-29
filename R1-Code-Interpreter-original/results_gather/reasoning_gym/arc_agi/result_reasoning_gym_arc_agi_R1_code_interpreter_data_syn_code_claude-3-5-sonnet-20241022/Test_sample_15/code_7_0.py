def transform_grid(input_grid):
    output_grid = [[2]*10 for _ in range(10)]
    
    # Copy original pattern
    for i in range(10):
        for j in range(10):
            output_grid[i][j] = input_grid[i][j]
    
    # Find pattern boundaries
    non_2_positions = []
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                non_2_positions.append((i, j))
    
    if not non_2_positions:
        return output_grid
    
    # Calculate center
    min_row = min(i for i, _ in non_2_positions)
    max_row = max(i for i, _ in non_2_positions)
    min_col = min(j for _, j in non_2_positions)
    max_col = max(j for _, j in non_2_positions)
    
    center_row = (min_row + max_row) // 2
    center_col = (min_col + max_col) // 2
    
    # Create reflections at specific distances
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                value = input_grid[i][j]
                
                # Calculate distances from center
                dist_row = i - center_row
                dist_col = j - center_col
                
                # Create reflection points at fixed distances
                reflection_points = []
                
                # If point is 2 units away from center horizontally
                if abs(dist_col) == 2:
                    reflection_points.append((i, center_col - dist_col))
                
                # If point is 2 units away from center vertically
                if abs(dist_row) == 2:
                    reflection_points.append((center_row - dist_row, j))
                    
                # If point is at both horizontal and vertical distance of 2
                if abs(dist_row) == 2 and abs(dist_col) == 2:
                    reflection_points.append((center_row - dist_row, center_col - dist_col))
                
                # Apply reflections
                for ref_row, ref_col in reflection_points:
                    if 0 <= ref_row < 10 and 0 <= ref_col < 10:
                        output_grid[ref_row][ref_col] = value
    
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