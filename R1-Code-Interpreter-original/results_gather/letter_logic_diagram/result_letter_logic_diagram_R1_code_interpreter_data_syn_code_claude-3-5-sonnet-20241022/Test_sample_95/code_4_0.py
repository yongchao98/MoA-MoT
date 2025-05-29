def print_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def get_minor_diagonal_value(grid):
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid_placement(grid, row, col, letter, initial_grid):
    # Must match pre-filled value
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
    if row + col == 6:  # Is on minor diagonal
        diag_value = get_minor_diagonal_value(grid)
        if diag_value and diag_value != letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, try only the diagonal value first
    if row + col == 6:
        diag_value = get_minor_diagonal_value(grid)
        if diag_value:
            if is_valid_placement(grid, row, col, diag_value, initial_grid):
                grid[row][col] = diag_value
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False

    # Try each possible letter
    letters = 'abcdefg'
    # If on minor diagonal without a fixed value, try one letter and stick with it
    if row + col == 6 and not diag_value:
        first_success = None
        for letter in letters:
            if is_valid_placement(grid, row, col, letter, initial_grid):
                first_success = letter
                break
        letters = first_success if first_success else letters

    for letter in letters:
        if is_valid_placement(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize the grid with pre-filled values
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
grid = [row[:] for row in initial_grid]

# Verify initial state and solve
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            if not is_valid_placement(grid, i, j, initial_grid[i][j], initial_grid):
                valid = False
                break
    if not valid:
        break

if not valid:
    print("Initial configuration is invalid!")
else:
    if solve(grid, initial_grid):
        print_grid(grid)
    else:
        print("No solution found")