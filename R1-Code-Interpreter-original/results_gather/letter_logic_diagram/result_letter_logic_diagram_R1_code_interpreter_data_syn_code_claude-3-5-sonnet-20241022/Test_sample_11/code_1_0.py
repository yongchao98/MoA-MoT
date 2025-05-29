def get_minor_diagonal_letter(grid):
    # Find the first non-empty cell in minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
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
            
    # Check minor diagonal
    if row + col == 6:  # Position is on minor diagonal
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def find_empty(grid):
    # Find first empty cell, prioritizing minor diagonal
    # First check minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Then check rest of grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, get the required letter
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter:
            letters = diag_letter
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
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

# Solve the puzzle
if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")