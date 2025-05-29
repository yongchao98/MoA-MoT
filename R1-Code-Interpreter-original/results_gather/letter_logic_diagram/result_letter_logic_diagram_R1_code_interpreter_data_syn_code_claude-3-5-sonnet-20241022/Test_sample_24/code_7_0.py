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

def is_valid(grid, r, c, val, initial):
    # Check initial constraints
    if initial[r][c] != '.' and initial[r][c] != val:
        return False

    # Check row
    if val in grid[r]:
        return False

    # Check column
    for i in range(7):
        if grid[i][c] == val:
            return False

    # Check minor diagonal requirement
    if r + c == 6:  # if on minor diagonal
        # Find any existing 'a' on the diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '.' and grid[i][j] != val:
                return False
            if initial[i][j] != '.' and initial[i][j] != val:
                return False

    return True

def find_next_empty(grid):
    # First try to fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '.':
            return (i, j)
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '.':
                return (i, j)
    return None

def solve(grid, initial):
    pos = find_next_empty(grid)
    if not pos:
        return True
    
    r, c = pos
    candidates = ['a'] if r + c == 6 else 'abcdefg'
    
    for val in candidates:
        if is_valid(grid, r, c, val, initial):
            grid[r][c] = val
            if solve(grid, initial):
                return True
            grid[r][c] = '.'
    
    return False

# Initialize grid
initial = get_initial_grid()
grid = [row[:] for row in initial]

# Solve
if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution exists")