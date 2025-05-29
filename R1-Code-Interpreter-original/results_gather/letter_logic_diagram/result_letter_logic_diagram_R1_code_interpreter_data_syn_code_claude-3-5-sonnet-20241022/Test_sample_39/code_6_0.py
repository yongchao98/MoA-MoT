def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid(grid, row, col, letter, initial):
    # First check: MUST match initial configuration
    if initial[row][col] != '' and initial[row][col] != letter:
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
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def solve(grid, initial, empty_cells):
    if not empty_cells:
        return True

    # Get next empty cell
    row, col = empty_cells[0]
    remaining_cells = empty_cells[1:]

    # If this is a pre-filled cell, we must use that value
    if initial[row][col] != '':
        if is_valid(grid, row, col, initial[row][col], initial):
            grid[row][col] = initial[row][col]
            if solve(grid, initial, remaining_cells):
                return True
        grid[row][col] = ''
        return False

    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, remaining_cells):
                return True
            grid[row][col] = ''

    return False

# Initial configuration
initial = [
    ['', 'b', '', '', '', '', ''],
    ['b', 'a', 'g', 'd', '', '', 'f'],
    ['', '', 'd', '', '', '', 'b'],
    ['g', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', ''],
    ['c', '', 'f', '', 'a', '', 'd'],
    ['e', 'f', '', 'a', 'g', '', '']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# Get list of all cells to fill
empty_cells = []
for i in range(7):
    for j in range(7):
        # First add pre-filled cells
        if initial[i][j] != '':
            empty_cells.append((i, j))
# Then add empty cells
for i in range(7):
    for j in range(7):
        if initial[i][j] == '':
            empty_cells.append((i, j))

# Try to solve
if solve(grid, initial, empty_cells):
    print_solution(grid)
else:
    print('<<<')
    print("No solution exists")
    print('>>>')