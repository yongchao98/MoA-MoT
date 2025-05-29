def print_grid(grid):
    for row in grid:
        print(row)

def is_valid_placement(grid, initial_grid, row, col, letter):
    # First check: must match pre-filled position if any
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Second check: minor diagonal must have same letter
    if row + col == 6:  # if on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    # Third check: row constraint
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Fourth check: column constraint
    for i in range(7):
        if grid[i][col] == letter:
            return False

    return True

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return i, j
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_sudoku(grid, initial_grid):
    # Find empty position
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos

    # If it's on the minor diagonal, get the required letter
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                if is_valid_placement(grid, initial_grid, row, col, grid[i][j]):
                    grid[row][col] = grid[i][j]
                    if solve_sudoku(grid, initial_grid):
                        return True
                    grid[row][col] = ''
                return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_placement(grid, initial_grid, row, col, letter):
            grid[row][col] = letter
            if solve_sudoku(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with pre-filled values
initial_grid = [
    ['', 'g', 'b', '', '', '', ''],
    ['', '', '', '', 'c', '', ''],
    ['b', '', '', '', '', 'a', ''],
    ['f', 'd', 'c', '', '', '', ''],
    ['', 'c', 'e', 'a', '', 'b', ''],
    ['', '', '', 'g', '', 'f', 'd'],
    ['e', 'a', 'g', '', 'f', '', '']
]

# Create working grid and copy initial values
grid = [row[:] for row in initial_grid]

# First, determine diagonal letter if possible
diag_letter = None
for i in range(7):
    j = 6 - i
    if initial_grid[i][j] != '':
        diag_letter = initial_grid[i][j]
        break

# Fill diagonal with known letter if found
if diag_letter:
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = diag_letter

if solve_sudoku(grid, initial_grid):
    print("<<<")
    for row in grid:
        print(",".join(row))
    print(">>>")
else:
    print("No solution exists")