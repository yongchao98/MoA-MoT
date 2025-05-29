def is_valid_partial(grid, row, col, num, constraints):
    # Check if num is already in the current row or column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, op) in constraints['horizontal']:
        if r == row:
            if c1 == col and grid[row][c2] != 0:
                if op == '<' and not (num < grid[row][c2]):
                    return False
                if op == '>' and not (num > grid[row][c2]):
                    return False
            if c2 == col and grid[row][c1] != 0:
                if op == '<' and not (grid[row][c1] < num):
                    return False
                if op == '>' and not (grid[row][c1] > num):
                    return False
    
    # Check vertical constraints
    for (c, r1, r2, op) in constraints['vertical']:
        if c == col:
            if r1 == row and grid[r2][col] != 0:
                if op == '∧' and not (num < grid[r2][col]):
                    return False
                if op == '∨' and not (num > grid[r2][col]):
                    return False
            if r2 == row and grid[r1][col] != 0:
                if op == '∧' and not (grid[r1][col] < num):
                    return False
                if op == '∨' and not (grid[r1][col] > num):
                    return False
    
    return True

def solve_futoshiki_backtrack(grid, constraints, row=0, col=0):
    if row == 7:
        return True
    
    next_row, next_col = (row, col + 1) if col < 6 else (row + 1, 0)
    
    if grid[row][col] != 0:
        return solve_futoshiki_backtrack(grid, constraints, next_row, next_col)
    
    for num in range(1, 8):
        if is_valid_partial(grid, row, col, num, constraints):
            grid[row][col] = num
            if solve_futoshiki_backtrack(grid, constraints, next_row, next_col):
                return True
            grid[row][col] = 0
    
    return False

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 5, 0, 0, 0, 0, 0],
        [4, 0, 0, 3, 0, 0, 0],
        [0, 3, 0, 5, 0, 1, 7],
        [0, 0, 0, 0, 0, 3, 0],
        [1, 0, 0, 0, 5, 0, 0],
        [0, 0, 0, 0, 0, 0, 4],
        [0, 0, 6, 2, 0, 0, 0]
    ]
    
    # Constraints
    constraints = {
        'horizontal': [
            (1, 0, 1, '<'), (1, 2, 3, '>'),
            (2, 3, 4, '<'), (2, 4, 5, '>'),
            (3, 3, 4, '<'),
            (4, 1, 2, '>'), (4, 4, 5, '>'),
            (5, 0, 1, '>'),
            (6, 1, 2, '<')
        ],
        'vertical': [
            (1, 0, 1, '∧'),
            (3, 1, 2, '∧'),
            (4, 2, 3, '∨'),
            (5, 3, 4, '∨'),
            (6, 4, 5, '∧')
        ]
    }
    
    if solve_futoshiki_backtrack(grid, constraints):
        for row in grid:
            print(row)
    else:
        print("No solution found")

solve_futoshiki()