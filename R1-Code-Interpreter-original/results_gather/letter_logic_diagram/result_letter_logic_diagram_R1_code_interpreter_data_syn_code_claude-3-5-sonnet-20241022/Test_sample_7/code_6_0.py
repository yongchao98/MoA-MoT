def get_minor_diagonal_letter(grid):
    # Find the first non-empty letter on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter, initial_grid):
    # Check 1: Must match initial grid position if specified
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check 2: Row uniqueness
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check 3: Column uniqueness
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # Check 4: Minor diagonal must have same letter
    if row + col == 6:  # if on minor diagonal
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
            return False

    return True

def solve(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # If this position is pre-filled in initial grid
    if initial_grid[row][col] != '':
        if not is_valid(grid, row, col, initial_grid[row][col], initial_grid):
            return False
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, row, col + 1)

    # If on minor diagonal
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter:
            if is_valid(grid, row, col, diag_letter, initial_grid):
                grid[row][col] = diag_letter
                if solve(grid, initial_grid, row, col + 1):
                    return True
            grid[row][col] = ''
            return False

    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''

    return False

# Initialize grids
initial_grid = [
    ['', 'g', 'f', '', 'e', '', ''],
    ['', 'f', 'a', '', 'c', '', 'b'],
    ['', '', '', 'c', 'd', '', ''],
    ['', '', 'c', 'd', 'b', '', ''],
    ['e', 'c', '', '', '', '', 'a'],
    ['', 'd', '', 'g', 'f', 'a', 'e'],
    ['', '', 'g', 'f', '', 'e', 'c']
]

grid = [['' for _ in range(7)] for _ in range(7)]

# First fill in all initial positions
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")