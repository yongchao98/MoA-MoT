def is_valid(grid, row, col, letter, initial_grid):
    # First check: must match initial grid positions
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check row uniqueness
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column uniqueness
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

def solve(grid, initial_grid):
    # First, verify and fill all initial positions
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                if grid[i][j] == '':
                    if not is_valid(grid, i, j, initial_grid[i][j], initial_grid):
                        return False
                    grid[i][j] = initial_grid[i][j]

    # Find empty position
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break

    if not empty:
        return True

    row, col = empty

    # If on minor diagonal, try to match existing diagonal letter first
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter:
            if is_valid(grid, row, col, diag_letter, initial_grid):
                grid[row][col] = diag_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False

    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
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

# Copy initial values
for i in range(7):
    for j in range(7):
        grid[i][j] = initial_grid[i][j]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")