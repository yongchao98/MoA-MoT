def find_non_2_region(grid):
    non_2_positions = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] != 2:
                non_2_positions.append((i, j))
    return non_2_positions

def transform_grid(input_grid):
    # Create output grid
    output_grid = [[2]*10 for _ in range(10)]
    
    # Copy original pattern
    for i in range(10):
        for j in range(10):
            output_grid[i][j] = input_grid[i][j]
    
    # Find the center of the pattern
    positions = find_non_2_region(input_grid)
    if not positions:
        return output_grid
        
    min_row = min(i for i, _ in positions)
    max_row = max(i for i, _ in positions)
    min_col = min(j for _, j in positions)
    max_col = max(j for _, j in positions)
    
    center_row = (min_row + max_row) // 2
    center_col = (min_col + max_col) // 2
    
    # Apply the transformation
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 2:
                # Calculate mirror positions
                dist_row = i - center_row
                dist_col = j - center_col
                
                # Mirror horizontally
                mirror_col = center_col - dist_col
                if 0 <= mirror_col < 10:
                    output_grid[i][mirror_col] = input_grid[i][j]
                
                # Mirror vertically
                mirror_row = center_row - dist_row
                if 0 <= mirror_row < 10:
                    output_grid[mirror_row][j] = input_grid[i][j]
                    
                # Mirror both (diagonal)
                if 0 <= mirror_row < 10 and 0 <= mirror_col < 10:
                    output_grid[mirror_row][mirror_col] = input_grid[i][j]

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

# Transform and print result
output_grid = [[2]*10 for _ in range(10)]
transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))