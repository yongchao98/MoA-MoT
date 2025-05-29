def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_diagonal_candidates(grid):
    # Find all letters that could be on the diagonal
    # by checking which letters are compatible with pre-filled positions
    candidates = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            candidates = {grid[i][6-i]}
            break
    return candidates

def is_valid(grid, row, col, letter, diag_letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check diagonal
    if row + col == 6 and letter != diag_letter:
        return False
    
    return True

def can_use_diagonal_letter(grid, diag_letter):
    # Check if this diagonal letter is compatible with pre-filled positions
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diag_letter:
            return False
    return True

def solve(grid, diag_letter, row=0, col=0):
    if row == 7:
        return True
    
    next_row = row + (col + 1) // 7
    next_col = (col + 1) % 7
    
    if grid[row][col] != '':
        if row + col == 6 and grid[row][col] != diag_letter:
            return False
        return solve(grid, diag_letter, next_row, next_col)
    
    letters = [diag_letter] if row + col == 6 else list('abcdefg')
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve(grid, diag_letter, next_row, next_col):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', '', '', '', '', 'a', ''],
    ['', 'e', '', 'd', '', '', ''],
    ['', 'c', '', 'a', '', 'g', 'f'],
    ['', '', '', 'b', '', 'f', 'e'],
    ['d', 'a', '', '', '', 'e', ''],
    ['', '', '', 'f', 'e', '', ''],
    ['b', 'g', 'f', 'e', '', '', '']
]

# Try each possible letter for the diagonal
for diag_letter in 'abcdefg':
    grid_copy = [row[:] for row in initial_grid]
    if can_use_diagonal_letter(grid_copy, diag_letter):
        if solve(grid_copy, diag_letter):
            print_grid(grid_copy)
            break