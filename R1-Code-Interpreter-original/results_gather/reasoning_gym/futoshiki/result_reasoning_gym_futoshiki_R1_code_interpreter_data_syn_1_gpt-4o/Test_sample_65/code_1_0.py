def is_valid_partial(grid, row, col, num, constraints):
    # Check if num is already in the current row or column
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, op) in constraints['horizontal']:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if op == '<' and not (num < grid[r][c2]):
                    return False
                if op == '>' and not (num > grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if op == '<' and not (grid[r][c1] < num):
                    return False
                if op == '>' and not (grid[r][c1] > num):
                    return False
    
    # Check vertical constraints
    for (c, r1, r2, op) in constraints['vertical']:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if op == '∧' and not (num < grid[r2][c]):
                    return False
                if op == '∨' and not (num > grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if op == '∧' and not (grid[r1][c] < num):
                    return False
                if op == '∨' and not (grid[r1][c] > num):
                    return False
    
    return True

def solve_futoshiki_backtrack(grid, constraints, row=0, col=0):
    if row == 5:
        return True
    if col == 5:
        return solve_futoshiki_backtrack(grid, constraints, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki_backtrack(grid, constraints, row, col + 1)
    
    for num in range(1, 6):
        if is_valid_partial(grid, row, col, num, constraints):
            grid[row][col] = num
            if solve_futoshiki_backtrack(grid, constraints, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

def solve_futoshiki():
    # Define the constraints
    constraints = {
        'horizontal': [
            (0, 0, 1, '<'), (0, 1, 2, '>'), (0, 2, 3, '<'), (0, 3, 4, '>'),
            (1, 0, 1, '<'), (1, 1, 2, '<'), (1, 2, 3, '>'),
            (2, 1, 2, '<'), (2, 2, 3, '<'),
            (3, 0, 1, '>'),
            (4, 0, 1, '>')
        ],
        'vertical': [
            (0, 0, 1, '∧'),
            (1, 1, 2, '∨'),
            (2, 2, 3, '∨'),
            (3, 3, 4, '∧')
        ]
    }
    
    # Initialize the grid with zeros
    grid = [[0] * 5 for _ in range(5)]
    
    # Solve the puzzle
    if solve_futoshiki_backtrack(grid, constraints):
        return grid
    else:
        return None

solution = solve_futoshiki()
for row in solution:
    print(row)