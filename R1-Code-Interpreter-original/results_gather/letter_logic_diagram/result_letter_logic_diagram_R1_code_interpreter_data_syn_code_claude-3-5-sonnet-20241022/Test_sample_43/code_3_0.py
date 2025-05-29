def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_complete_valid(grid, initial):
    # Check all initial values
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    
    # Check rows and columns
    for i in range(7):
        row_vals = set(grid[i])
        col_vals = set(grid[r][i] for r in range(7))
        if len(row_vals) != 7 or len(col_vals) != 7:
            return False
    
    # Check minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1 or diagonal[0] != 'e':
        return False
    
    return True

def get_next_empty(grid, initial, r, c):
    while r < 7:
        while c < 7:
            if initial[r][c] == '':
                return r, c
            c += 1
        r += 1
        c = 0
    return None

def is_valid_move(grid, r, c, val):
    # Check row
    if val in grid[r]:
        return False
    
    # Check column
    if val in [grid[i][c] for i in range(7)]:
        return False
    
    # Check minor diagonal
    if r + c == 6:  # If on minor diagonal
        if val != 'e':
            return False
    
    return True

def solve(grid, initial):
    # Find next empty cell
    pos = get_next_empty(grid, initial, 0, 0)
    if not pos:
        return is_complete_valid(grid, initial)
    
    r, c = pos
    
    # If on minor diagonal, only try 'e'
    if r + c == 6:
        if is_valid_move(grid, r, c, 'e'):
            grid[r][c] = 'e'
            if solve(grid, initial):
                return True
        grid[r][c] = ''
        return False
    
    # Try each possible value
    for val in 'abcdefg':
        if is_valid_move(grid, r, c, val):
            grid[r][c] = val
            if solve(grid, initial):
                return True
            grid[r][c] = ''
    
    return False

# Initialize grid with initial values
initial = [
    ['d','a','','','g','',''],
    ['','','','g','e','',''],
    ['f','','g','','c','',''],
    ['','g','','','','a',''],
    ['g','e','','','','f',''],
    ['e','','','','','b',''],
    ['','d','','f','b','','e']
]

# Create working grid and copy initial values
grid = [['' for _ in range(7)] for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            grid[i][j] = initial[i][j]

# Solve
if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution found")