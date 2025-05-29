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

def verify_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '.' and grid[i][j] != initial[i][j]:
                return False
    return True

def is_valid(grid, r, c, val, initial):
    # Check initial constraints
    if initial[r][c] != '.' and initial[r][c] != val:
        return False

    # Check row
    for j in range(7):
        if grid[r][j] == val:
            return False

    # Check column
    for i in range(7):
        if grid[i][c] == val:
            return False

    # Check minor diagonal
    if r + c == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '.' and grid[i][j] != val:
                return False
            if initial[i][j] != '.' and initial[i][j] != val:
                return False

    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '.':
                return i, j
    return None

def solve(grid, initial):
    empty = find_empty(grid)
    if not empty:
        return verify_initial(grid, initial)
    
    r, c = empty
    
    # If on minor diagonal, find the required value
    if r + c == 6:
        diag_val = None
        for i in range(7):
            if initial[i][6-i] != '.':
                diag_val = initial[i][6-i]
                break
        if diag_val:
            if is_valid(grid, r, c, diag_val, initial):
                grid[r][c] = diag_val
                if solve(grid, initial):
                    return True
                grid[r][c] = '.'
            return False

    for val in 'abcdefg':
        if is_valid(grid, r, c, val, initial):
            grid[r][c] = val
            if solve(grid, initial):
                return True
            grid[r][c] = '.'
    
    return False

# Initialize grids
initial_grid = get_initial_grid()
working_grid = [row[:] for row in initial_grid]

# Solve
if solve(working_grid, initial_grid):
    print_grid(working_grid)
else:
    print("No solution exists")