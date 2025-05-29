def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_minor_diagonal_cells():
    return [(i, 6-i) for i in range(7)]

def get_required_diagonal_letter(grid):
    diag_cells = get_minor_diagonal_cells()
    letters = {}
    # Check pre-filled diagonal cells
    for i, j in diag_cells:
        if grid[i][j] != '':
            return grid[i][j]
    return None

def is_valid(grid, row, col, letter, diag_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal requirement
    if row + col == 6 and letter != diag_letter:
        return False
    
    return True

def find_empty(grid):
    # First check non-diagonal cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                return (i, j)
    
    # Then check diagonal cells
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    return None

def solve_grid(grid, diag_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If this is a diagonal position, only try the diagonal letter
    if row + col == 6:
        letters = diag_letter
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve_grid(grid, diag_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', 'd', 'g', '', 'b', 'f'],
    ['', 'd', '', '', '', 'f', ''],
    ['d', '', 'g', '', '', 'f', ''],
    ['g', '', 'b', 'f', 'e', 'c', 'd'],
    ['a', 'b', 'f', '', 'c', 'd', ''],
    ['', '', '', 'c', 'd', 'g', ''],
    ['f', '', 'c', 'd', 'g', 'a', 'b']
]

# First, determine what letter must be on the minor diagonal
diag_letter = get_required_diagonal_letter(initial_grid)
if not diag_letter:
    # If no diagonal letter is pre-filled, try each possibility
    for letter in 'abcdefg':
        grid_copy = [row[:] for row in initial_grid]
        # Fill diagonal with current letter
        for i, j in get_minor_diagonal_cells():
            if grid_copy[i][j] == '':
                if is_valid(grid_copy, i, j, letter, letter):
                    grid_copy[i][j] = letter
                else:
                    break
        else:
            if solve_grid(grid_copy, letter):
                print_grid(grid_copy)
                exit()
else:
    if solve_grid(initial_grid, diag_letter):
        print_grid(initial_grid)
    else:
        print("No solution found")