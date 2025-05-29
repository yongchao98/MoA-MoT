def is_valid_move(grid, row, col, num, constraints):
    # Check if num is not in the current row and column
    if num in grid[row] or num in [grid[r][col] for r in range(6)]:
        return False
    
    # Check horizontal constraints
    for (r, c1, c2) in constraints['horizontal']:
        if r == row:
            if c1 == col and grid[row][c2] != 0 and not (num > grid[row][c2]):
                return False
            if c2 == col and grid[row][c1] != 0 and not (grid[row][c1] > num):
                return False
    
    # Check vertical constraints
    for (c, r1, r2) in constraints['vertical']:
        if c == col:
            if r1 == row and grid[r2][col] != 0 and not (num < grid[r2][col]):
                return False
            if r2 == row and grid[r1][col] != 0 and not (grid[r1][col] < num):
                return False
    
    return True

def solve_futoshiki(grid, constraints, row=0, col=0):
    if row == 6:
        return True
    
    if col == 6:
        return solve_futoshiki(grid, constraints, row + 1, 0)
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, constraints, row, col + 1)
    
    for num in range(1, 7):
        if is_valid_move(grid, row, col, num, constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, constraints, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

def futoshiki_solver():
    grid = [
        [0, 3, 0, 6, 5, 0],
        [0, 0, 0, 2, 0, 6],
        [0, 0, 4, 0, 6, 2],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 5],
        [0, 5, 6, 3, 0, 0]
    ]
    
    constraints = {
        'horizontal': [(2, 4, 5), (4, 0, 1)],
        'vertical': [(5, 2, 3), (5, 3, 5)]
    }
    
    if solve_futoshiki(grid, constraints):
        for row in grid:
            print(row)
    else:
        print("No solution found")

futoshiki_solver()