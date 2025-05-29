def is_valid(grid, row, col, num, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7) if grid[i][col] != 0]:
        return False
    
    # Check vertical constraints
    for (r1, c1, r2, c2) in v_constraints:
        if row == r1 and col == c1 and grid[r2][c2] != 0:
            if num >= grid[r2][c2]:
                return False
        if row == r2 and col == c2 and grid[r1][c1] != 0:
            if grid[r1][c1] >= num:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 8):
        if is_valid(grid, row, col, num, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 0, 0, 7, 0, 0],
    [3, 2, 0, 0, 0, 0, 0],
    [0, 0, 2, 5, 0, 0, 7],
    [0, 0, 6, 0, 1, 5, 0],
    [4, 0, 0, 2, 6, 1, 0],
    [7, 0, 0, 0, 0, 3, 0],
    [0, 5, 0, 0, 0, 2, 4]
]

# Vertical constraints (row1, col1, row2, col2) where row1,col1 < row2,col2
v_constraints = [
    (1, 4, 2, 4),  # Second constraint
    (3, 4, 4, 4)   # Fourth constraint
]

if solve_futoshiki(grid, v_constraints):
    for row in grid:
        print(' '.join(map(str, row)))
else:
    print("No solution exists")