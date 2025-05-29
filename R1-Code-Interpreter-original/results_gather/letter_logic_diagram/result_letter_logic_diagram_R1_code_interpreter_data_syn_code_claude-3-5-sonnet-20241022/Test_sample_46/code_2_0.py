def get_minor_diagonal_letter(grid):
    # Get the letter that appears on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # If this is on minor diagonal
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def find_empty(grid):
    # First priority: fill minor diagonal if not complete
    diag_letter = get_minor_diagonal_letter(grid)
    if diag_letter:
        for i in range(7):
            if grid[i][6-i] == '':
                return i, 6-i
    
    # Second priority: fill remaining cells
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
    
    # If on minor diagonal, only use the diagonal letter if known
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        letters = [diag_letter] if diag_letter else 'abcdefg'
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'g', 'e', '', 'b', '', ''],
    ['g', 'e', 'f', 'b', '', '', ''],
    ['', '', '', '', '', '', 'g'],
    ['f', '', 'd', '', '', 'g', 'e'],
    ['', 'd', '', '', '', '', ''],
    ['d', 'c', '', '', '', '', 'b'],
    ['c', '', 'g', '', '', 'b', '']
]

# First, try each possible letter for minor diagonal
for diag_letter in 'abcdefg':
    # Create a copy of initial grid
    test_grid = [row[:] for row in initial_grid]
    
    # Check if this diagonal letter is possible
    valid = True
    for i in range(7):
        if test_grid[i][6-i] != '' and test_grid[i][6-i] != diag_letter:
            valid = False
            break
    
    if not valid:
        continue
        
    # Fill empty diagonal positions
    for i in range(7):
        if test_grid[i][6-i] == '':
            test_grid[i][6-i] = diag_letter
    
    # Try to solve with this diagonal
    if solve_grid(test_grid):
        for row in test_grid:
            print(','.join(row))
        break