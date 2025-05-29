def print_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid(grid, row, col, letter, initial_grid):
    # Check if we're trying to change a pre-filled cell
    if initial_grid[row][col] != '' and letter != initial_grid[row][col]:
        return False
    
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def solve_grid(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # If this is a pre-filled cell, skip to next cell
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve_grid(grid, initial_grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with pre-filled values
initial_grid = [
    ['', '', '', '', 'g', '', ''],
    ['', '', '', '', '', '', ''],
    ['e', '', 'g', 'c', '', '', 'b'],
    ['', 'g', '', 'd', 'f', '', 'e'],
    ['g', 'c', '', '', 'b', '', 'a'],
    ['c', 'd', '', '', 'e', '', 'g'],
    ['', 'f', '', 'e', '', '', '']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# First, verify all pre-filled cells
valid_initial = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            if not is_valid(grid, i, j, initial_grid[i][j], initial_grid):
                valid_initial = False
                break
            grid[i][j] = initial_grid[i][j]
    if not valid_initial:
        break

if valid_initial and solve_grid(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution found")