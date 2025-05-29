def analyze_grid_features(input_grid):
    rows = input_grid.split('\n')
    grid = [list(map(int, row.split())) for row in rows]
    height = len(grid)
    width = len(grid[0])
    
    # Count 2's in different regions
    twos_center = 0
    twos_edges = 0
    
    # Define center region
    center_start_h = height // 3
    center_end_h = 2 * height // 3
    center_start_w = width // 3
    center_end_w = 2 * width // 3
    
    for i in range(height):
        for j in range(width):
            if grid[i][j] == 2:
                if (center_start_h <= i <= center_end_h and 
                    center_start_w <= j <= center_end_w):
                    twos_center += 1
                else:
                    twos_edges += 1
    
    return {
        'twos_center': twos_center,
        'twos_edges': twos_edges,
        'total_twos': twos_center + twos_edges
    }

test_input = """0 0 3 2 3 3 3 3
0 0 3 3 3 3 3 3
3 3 3 3 3 2 2 3
3 2 2 3 3 2 2 3
3 2 2 3 0 3 3 3"""

features = analyze_grid_features(test_input)
print("Test input features:", features)