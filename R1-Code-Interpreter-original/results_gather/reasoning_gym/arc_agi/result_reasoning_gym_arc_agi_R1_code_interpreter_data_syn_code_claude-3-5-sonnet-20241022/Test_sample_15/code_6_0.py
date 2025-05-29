def transform_grid(input_grid):
    output_grid = [[2]*10 for _ in range(10)]
    
    # Copy original pattern
    for i in range(10):
        for j in range(10):
            output_grid[i][j] = input_grid[i][j]
    
    # Find the center of the pattern
    center_row = 4  # Based on all examples, the pattern is centered around row 4
    center_col = 5  # and column 5
    
    # Create reflections at specific positions
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                # Calculate distance from center
                dist_row = i - center_row
                dist_col = j - center_col
                
                # Create reflection points
                reflection_points = []
                
                # If point is on one side of center, reflect to other side
                if dist_col != 0:  # Horizontal reflection
                    new_col = center_col - dist_col
                    if 0 <= new_col < 10:
                        reflection_points.append((i, new_col))
                
                if dist_row != 0:  # Vertical reflection
                    new_row = center_row - dist_row
                    if 0 <= new_row < 10:
                        reflection_points.append((new_row, j))
                        # Also add diagonal reflection
                        if dist_col != 0:
                            new_col = center_col - dist_col
                            if 0 <= new_col < 10:
                                reflection_points.append((new_row, new_col))
                
                # Apply reflections
                value = input_grid[i][j]
                for ref_row, ref_col in reflection_points:
                    if input_grid[ref_row][ref_col] == 2:  # Only replace if it's empty (2)
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