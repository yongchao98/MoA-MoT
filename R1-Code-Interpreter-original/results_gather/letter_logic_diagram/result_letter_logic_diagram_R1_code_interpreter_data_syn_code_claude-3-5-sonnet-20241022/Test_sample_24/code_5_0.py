def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_initial_grid():
    return [
        ['f', 'd', 'a', '.', '.', 'g', '.'],
        ['d', '.', '.', 'b', '.', '.', 'f'],
        ['.', 'c', 'b', '.', 'e', '.', '.'],
        ['.', '.', 'g', 'e', '.', '.', '.'],
        ['.', 'g', 'e', 'f', '.', 'a', 'c'],
        ['.', 'e', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'a', '.', '.', '.']
    ]

def check_row_col(grid, r, c, val):
    # Check row
    if val in grid[r]:
        return False
    # Check column
    if val in [grid[i][c] for i in range(7)]:
        return False
    return True

def get_diagonal_value(grid, initial_grid):
    # Find the first non-empty value on the minor diagonal
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] != '.':
            return initial_grid[i][j]
    return None

def is_valid_move(grid, initial_grid, r, c, val):
    # Check initial constraints
    if initial_grid[r][c] != '.' and initial_grid[r][c] != val:
        return False
    
    # Check row and column
    if not check_row_col(grid, r, c, val):
        return False
    
    # Check minor diagonal
    if r + c == 6:  # If on minor diagonal
        diag_val = get_diagonal_value(grid, initial_grid)
        if diag_val and val != diag_val:
            return False
    
    return True

def solve_puzzle():
    initial_grid = get_initial_grid()
    grid = [['.' for _ in range(7)] for _ in range(7)]
    
    # First, fill in all known values from initial grid
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '.':
                grid[i][j] = initial_grid[i][j]
    
    # Find diagonal value
    diag_val = get_diagonal_value(grid, initial_grid)
    if not diag_val:
        diag_val = 'a'  # If no diagonal value is given, use 'a'
    
    # Fill diagonal first
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '.':
            grid[i][j] = diag_val
    
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
        
        r, c = pos // 7, pos % 7
        
        # Skip if cell is already filled
        if grid[r][c] != '.':
            return backtrack(pos + 1)
        
        # Try each possible value
        for val in 'abcdefg':
            if is_valid_move(grid, initial_grid, r, c, val):
                grid[r][c] = val
                if backtrack(pos + 1):
                    return True
                grid[r][c] = '.'
        
        return False
    
    if backtrack():
        return grid
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    print_grid(solution)
else:
    print("No solution exists")