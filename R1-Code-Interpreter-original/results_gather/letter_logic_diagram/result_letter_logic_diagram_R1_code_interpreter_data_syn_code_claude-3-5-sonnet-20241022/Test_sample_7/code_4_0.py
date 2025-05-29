def verify_initial_positions(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def is_valid(grid, row, col, letter, initial_grid):
    # First priority: must match initial grid
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

    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def find_empty(grid, initial_grid):
    # First fill positions that are in initial grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and initial_grid[i][j] != '':
                return (i, j, True)  # True indicates this is an initial position
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)
    return None

def solve(grid, initial_grid):
    empty = find_empty(grid, initial_grid)
    if not empty:
        return verify_initial_positions(grid, initial_grid)

    row, col, is_initial = empty

    if is_initial:
        # For initial positions, only try the specified letter
        letter = initial_grid[row][col]
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    else:
        # For other positions
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

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")