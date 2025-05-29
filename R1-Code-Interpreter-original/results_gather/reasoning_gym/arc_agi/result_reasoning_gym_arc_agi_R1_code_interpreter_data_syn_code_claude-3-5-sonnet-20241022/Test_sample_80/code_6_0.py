def analyze_cluster_distribution(grid_str):
    rows = grid_str.split('\n')
    grid = [list(map(int, row.split())) for row in rows]
    height = len(grid)
    width = len(grid[0])
    
    # Analyze horizontal distribution of 2's
    left_half = sum(1 for i in range(height) 
                   for j in range(width//2) if grid[i][j] == 2)
    right_half = sum(1 for i in range(height) 
                    for j in range(width//2, width) if grid[i][j] == 2)
    
    # Analyze vertical distribution
    top_half = sum(1 for i in range(height//2) 
                  for j in range(width) if grid[i][j] == 2)
    bottom_half = sum(1 for i in range(height//2, height) 
                     for j in range(width) if grid[i][j] == 2)
    
    return {
        'left_right_ratio': f"{left_half}:{right_half}",
        'top_bottom_ratio': f"{top_half}:{bottom_half}",
        'total_twos': left_half + right_half
    }

test_input = """0 0 3 2 3 3 3 3
0 0 3 3 3 3 3 3
3 3 3 3 3 2 2 3
3 2 2 3 3 2 2 3
3 2 2 3 0 3 3 3"""

distribution = analyze_cluster_distribution(test_input)
print("Test input distribution:", distribution)