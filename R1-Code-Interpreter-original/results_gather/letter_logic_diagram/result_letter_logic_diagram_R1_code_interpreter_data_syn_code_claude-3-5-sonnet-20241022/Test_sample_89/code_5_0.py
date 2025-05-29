def print_grid(grid):
    print("<<<")
    for row in grid:
        print(','.join(row))
    print(">>>")

initial_grid = [
    ['', '', 'a', 'g', '', 'b', 'c'],
    ['', 'a', '', '', '', 'c', 'd'],
    ['a', 'g', 'e', 'b', '', '', ''],
    ['g', '', 'b', '', 'd', 'f', 'a'],
    ['e', 'b', '', '', 'f', '', 'g'],
    ['', '', 'd', 'f', 'a', '', ''],
    ['', 'd', 'f', 'a', 'g', 'e', 'b']
]

def is_valid(grid, row, col, letter, initial_grid):
    # Check if this conflicts with initial grid
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
            if initial_grid[i][j] != '' and initial_grid[i][j] != letter:
                return False

    return True

def find_minor_diagonal_letter(initial_grid):
    for i in range(7):
        if initial_grid[i][6-i] != '':
            return initial_grid[i][6-i]
    return 'f'  # default if none found

def solve(grid, initial_grid):
    # Find empty cell
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
    
    # If on minor diagonal, use only the required letter
    if row + col == 6:
        letter = find_minor_diagonal_letter(initial_grid)
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
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

# Create working grid and copy initial values
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

# First fill minor diagonal with the same letter
diag_letter = find_minor_diagonal_letter(initial_grid)
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = diag_letter

if solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")