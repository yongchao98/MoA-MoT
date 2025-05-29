def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_grid_valid(grid):
    # Check rows
    for row in grid:
        if row.count('') > 0:
            return False
        if len(set(row)) != 7:
            return False
    
    # Check columns
    for col in range(7):
        column = [grid[row][col] for row in range(7)]
        if len(set(column)) != 7:
            return False
    
    # Check diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def is_valid_move(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # Check diagonal constraint
    if row + col == 6:
        # If there's already a letter on diagonal, new letter must match
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return check_grid_valid(grid)
    
    row, col = empty
    
    # If this is a diagonal position, we must use the diagonal letter
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter:
            if is_valid_move(grid, row, col, diag_letter):
                grid[row][col] = diag_letter
                if solve(grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
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

# First, determine what letter must be on the diagonal
diag_letter = None
for i in range(7):
    if initial_grid[i][6-i] != '':
        diag_letter = initial_grid[i][6-i]
        break

if diag_letter:
    # Pre-fill all empty diagonal positions
    for i in range(7):
        if initial_grid[i][6-i] == '':
            initial_grid[i][6-i] = diag_letter

if solve(initial_grid):
    print_grid(initial_grid)
else:
    print("No solution found")