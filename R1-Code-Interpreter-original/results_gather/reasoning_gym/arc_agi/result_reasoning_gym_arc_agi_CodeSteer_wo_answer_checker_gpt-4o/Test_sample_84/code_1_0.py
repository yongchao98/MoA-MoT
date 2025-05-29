def count_nines_in_grid(grid):
    n = len(grid)
    m = len(grid[0]) if n > 0 else 0
    
    # Initialize counts
    max_nines = 0
    
    # Check rows
    for row in grid:
        max_nines = max(max_nines, row.count(9))
    
    # Check columns
    for col in range(m):
        col_count = sum(1 for row in grid if row[col] == 9)
        max_nines = max(max_nines, col_count)
    
    # Check main diagonal
    main_diag_count = sum(1 for i in range(min(n, m)) if grid[i][i] == 9)
    max_nines = max(max_nines, main_diag_count)
    
    # Check anti-diagonal
    anti_diag_count = sum(1 for i in range(min(n, m)) if grid[i][m - i - 1] == 9)
    max_nines = max(max_nines, anti_diag_count)
    
    return max_nines

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 6, 6, 6, 9, 9, 9],
    [9, 9, 9, 6, 9, 6, 6, 6, 9],
    [9, 6, 6, 6, 9, 9, 9, 6, 6],
    [6, 6, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 6, 6, 6, 9, 9],
    [6, 6, 9, 9, 6, 9, 6, 6, 6],
    [9, 6, 6, 6, 6, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9]
]

output_length = count_nines_in_grid(input_grid)
output_grid = [9] * output_length
print("<<<" + " ".join(map(str, output_grid)) + ">>>")