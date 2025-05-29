def is_valid_move(grid, row, col, num, constraints):
    # Check if num is not in the current row and column
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check constraints
    for (r1, c1, r2, c2, op) in constraints:
        if (r1 == row and c1 == col) or (r2 == row and c2 == col):
            val1 = grid[r1][c1] if (r1 != row or c1 != col) else num
            val2 = grid[r2][c2] if (r2 != row or c2 != col) else num
            if op == '<' and not (val1 < val2):
                return False
            if op == '>' and not (val1 > val2):
                return False
    return True

def solve_futoshiki_backtrack(grid, constraints, row=0, col=0):
    if row == 6:
        return True
    if col == 6:
        return solve_futoshiki_backtrack(grid, constraints, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki_backtrack(grid, constraints, row, col + 1)
    
    for num in range(1, 7):
        if is_valid_move(grid, row, col, num, constraints):
            grid[row][col] = num
            if solve_futoshiki_backtrack(grid, constraints, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

def solve_futoshiki():
    # Initial grid with given numbers
    grid = [
        [0, 0, 0, 0, 6, 2],
        [6, 0, 2, 0, 0, 4],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [5, 0, 0, 0, 3, 6],
        [0, 5, 0, 3, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 1, 1, 1, '>'),
        (0, 4, 1, 4, '>'),
        (2, 5, 3, 5, '>'),
        (4, 2, 4, 3, '<'),
        (5, 0, 5, 1, '<'),
        (5, 1, 5, 2, '>')
    ]
    
    if solve_futoshiki_backtrack(grid, constraints):
        return grid
    else:
        return None

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")