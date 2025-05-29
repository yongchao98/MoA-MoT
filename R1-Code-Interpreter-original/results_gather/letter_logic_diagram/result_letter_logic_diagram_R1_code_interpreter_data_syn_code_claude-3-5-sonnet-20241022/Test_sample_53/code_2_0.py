def get_minor_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter):
    # If this is a minor diagonal position, it must match the minor diagonal letter
    if row + col == 6:
        minor_letter = get_minor_diagonal_letter(grid)
        if minor_letter and letter != minor_letter:
            return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    # First priority: fill minor diagonal
    minor_letter = get_minor_diagonal_letter(grid)
    if minor_letter:
        for i in range(7):
            if grid[i][6-i] == '':
                return i, 6-i
    
    # Second priority: fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If we're on minor diagonal and we already have a letter there
    if row + col == 6:
        minor_letter = get_minor_diagonal_letter(grid)
        if minor_letter:
            if is_valid(grid, row, col, minor_letter):
                grid[row][col] = minor_letter
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
initial_grid = [
    ['', '', 'f', '', '', '', ''],
    ['', '', 'd', 'e', 'g', '', ''],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', '', 'b', '', 'f', ''],
    ['', 'c', '', 'a', 'f', '', ''],
    ['', '', '', 'f', 'd', 'e', 'g']
]

# First, determine what letter should be on minor diagonal
minor_letter = get_minor_diagonal_letter(grid)
if not minor_letter:
    # If no letter is preset on minor diagonal, try each possibility
    for letter in 'abcdefg':
        grid[0][6] = letter  # Set first minor diagonal position
        if solve_grid(grid):
            break
        grid[0][6] = ''
else:
    solve_grid(grid)

for row in initial_grid:
    print(','.join(row))