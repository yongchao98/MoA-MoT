def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # If on minor diagonal, must match diagonal_letter
    if row + col == 6 and diagonal_letter is not None:
        if letter != diagonal_letter:
            return False
            
    return True

def find_diagonal_letter(grid):
    # Check pre-filled cells on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then fill rest of grid
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
    letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['f','','e','','','',''],
    ['','e','','d','','','f'],
    ['','','','','g','f','a'],
    ['','d','b','','','a',''],
    ['d','b','','','a','e',''],
    ['b','g','','a','e','c','d'],
    ['','f','a','','','d','']
]

# First, find the diagonal letter
diagonal_letter = find_diagonal_letter(grid)
if diagonal_letter is None:
    diagonal_letter = 'c'  # We'll try 'c' as it appears in one pre-filled position

# Fill all minor diagonal positions with the diagonal letter
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = diagonal_letter

if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")