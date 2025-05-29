def print_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def get_minor_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter, initial_grid):
    # Must match pre-filled cell if exists
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # Check minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
            return False
        
    return True

def solve_grid(grid, initial_grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
        
    row = pos // 7
    col = pos % 7
    
    # If pre-filled, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve_grid(grid, initial_grid, pos + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid, pos + 1):
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

# First, copy all pre-filled values and verify they don't violate constraints
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

if not valid_initial:
    print("Initial configuration is invalid!")
else:
    if solve_grid(grid, initial_grid):
        print_grid(grid)
    else:
        print("No solution found")