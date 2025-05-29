def find_pattern_center(grid):
    non_2_positions = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] != 2:
                non_2_positions.append((i, j))
    
    # Find center of non-2 positions
    if non_2_positions:
        min_x = min(x for x, _ in non_2_positions)
        max_x = max(x for x, _ in non_2_positions)
        min_y = min(y for _, y in non_2_positions)
        max_y = max(y for _, y in non_2_positions)
        center_x = (min_x + max_x) // 2
        center_y = (min_y + max_y) // 2
        return center_x, center_y
    return None

def transform_grid(input_grid):
    output_grid = [[2]*10 for _ in range(10)]
    center_x, center_y = find_pattern_center(input_grid)
    
    # Copy the original pattern
    for i in range(len(input_grid)):
        for j in range(len(input_grid[0])):
            if input_grid[i][j] != 2:
                output_grid[i][j] = input_grid[i][j]
    
    # Apply reflective transformations
    for i in range(len(input_grid)):
        for j in range(len(input_grid[0])):
            if input_grid[i][j] != 2:
                # Reflect horizontally
                mirror_j = 2*center_y - j
                if 0 <= mirror_j < 10:
                    output_grid[i][mirror_j] = input_grid[i][j]
                
                # Reflect vertically
                mirror_i = 2*center_x - i
                if 0 <= mirror_i < 10:
                    output_grid[mirror_i][j] = input_grid[i][j]
                    
                # Reflect diagonally
                if 0 <= mirror_i < 10 and 0 <= mirror_j < 10:
                    output_grid[mirror_i][mirror_j] = input_grid[i][j]

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

output_grid = [[2]*10 for _ in range(10)]
transform_grid(input_grid)

# Print the result
for row in output_grid:
    print(' '.join(map(str, row)))