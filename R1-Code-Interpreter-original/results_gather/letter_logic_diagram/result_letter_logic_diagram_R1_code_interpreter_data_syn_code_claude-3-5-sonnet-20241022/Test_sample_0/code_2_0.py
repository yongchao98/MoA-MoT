def print_grid(grid):
    for row in grid:
        print(','.join(str(x) if x else '' for x in row))

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal
    if row + col == 6:  # Position is on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j, True)  # True indicates diagonal position
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col, is_diagonal = empty
    
    # If filling a diagonal position, first try the letter that's already on the diagonal
    if is_diagonal:
        diag_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diag_letter = grid[i][j]
                break
        if diag_letter:
            if is_valid(grid, row, col, diag_letter):
                grid[row][col] = diag_letter
                if solve_grid(grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', 'e', '', '', '', ''],
    ['g', '', 'e', 'f', '', '', 'c'],
    ['', '', '', 'a', '', '', ''],
    ['', '', '', '', '', '', 'e'],
    ['', 'a', 'd', 'c', '', '', ''],
    ['', 'd', 'c', 'g', '', '', 'b'],
    ['d', 'c', 'g', '', 'f', 'b', '']
]

# First, determine what letter must be on the minor diagonal (must be 'e' due to row 4, col 6)
diagonal_letter = 'e'

# Fill all diagonal positions with 'e'
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = diagonal_letter

if solve_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(x for x in row))
    print('>>>')
else:
    print("No solution exists")