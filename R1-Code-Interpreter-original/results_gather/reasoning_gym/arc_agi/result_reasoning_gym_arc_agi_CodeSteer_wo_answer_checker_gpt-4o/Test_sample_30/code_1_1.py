def analyze_grid(input_grid):
    # Analyze the grid for patterns
    # For simplicity, let's check for symmetry and repeating sequences
    n = len(input_grid)
    m = len(input_grid[0])
    
    # Check for horizontal symmetry
    horizontal_symmetry = all(input_grid[i] == input_grid[n-i-1] for i in range(n//2))
    
    # Check for vertical symmetry
    vertical_symmetry = all(row[j] == row[m-j-1] for row in input_grid for j in range(m//2))
    
    # Check for repeating sequences in rows
    repeating_rows = any(row.count(row[0]) == len(row) for row in input_grid)
    
    # Check for repeating sequences in columns
    repeating_columns = any(all(input_grid[i][j] == input_grid[0][j] for i in range(n)) for j in range(m))
    
    # Print findings
    print("Horizontal Symmetry:", horizontal_symmetry)
    print("Vertical Symmetry:", vertical_symmetry)
    print("Repeating Rows:", repeating_rows)
    print("Repeating Columns:", repeating_columns)

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 1, 9, 9, 9, 0, 0, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 1, 9, 1, 9, 9, 0, 9, 9, 9],
    [9, 9, 9, 8, 9, 8, 9, 9, 1, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 8, 9, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 7, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 5, 6, 6, 6, 8, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 5, 9, 5, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 5, 5, 5, 9, 9, 9, 9],
    [6, 7, 6, 6, 6, 0, 6, 9, 9, 9, 5, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9]
]

analyze_grid(input_grid)